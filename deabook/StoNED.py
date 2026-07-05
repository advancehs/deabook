# import dependencies
import numpy as np
import pandas
import pandas as pd
import math
import scipy.stats as stats
import scipy.optimize as opt
from dataclasses import dataclass
from .utils import tools
from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RED_MOM, RED_QLE, RED_KDE
from pyomo.environ import   exp
from math import sqrt, pi, log
from scipy.signal import convolve, correlate


@dataclass
class StoNEDContext:
    """Small adapter around model variants used by StoNED."""

    model: object
    family: str
    basexy: object = None
    gb: object = None
    gddf_er: object = None
    gddf: object = None
    gresidual: object = None
    basexy_old: object = None

    @classmethod
    def from_model(cls, model):
        name = model.__class__.__name__
        if not hasattr(model, "get_residual"):
            raise ValueError(f"Unsupported StoNED model type without residuals: {name}")
        if not hasattr(model, "fun"):
            raise ValueError(f"Unsupported StoNED model type without production/cost flag: {name}")

        return cls(
            model=model,
            family=name,
            basexy=getattr(model, "basexyb", getattr(model, "basexy", None)),
            gb=getattr(model, "gb", None),
            gddf_er=getattr(model, "gddf_er", None),
            gddf=getattr(model, "gddf", None),
            gresidual=getattr(model, "gresidual", None),
            basexy_old=getattr(model, "basexyb_old", None),
        )

    @property
    def is_ddf(self):
        return "DDF" in self.family or self.family in ["CNLSDDF", "CNLSDDFG", "CNLSDDFweak", "CNLSDDFweakG"]

    @property
    def is_meta_ddf(self):
        return self.family in ["CNLSSDweakmeta", "CNLSDDFweakmeta"]


class StoNED:
    """Stochastic nonparametric envelopment of data (StoNED)
    """

    def __init__(self, model):
        """StoNED
        model: The input model for residual decomposition
        """

        # print(model.__class__.__name__) # CNLSSD


        self.model = model
        self.context = StoNEDContext.from_model(model)
        self.basexy = self.context.basexy
        self.gb = self.context.gb
        self.gddf_er = self.context.gddf_er
        self.gddf = self.context.gddf
        self.gresidual = self.context.gresidual
        self.basexy_old = self.context.basexy_old

        self.y = getattr(model, "y", None)
        # self.data = model.data
        # self.sent = model.sent
        # self.z = model.z
        self.gy = getattr(model, "gy", None)
        self.gx = getattr(model, "gx", None)
        # self.rts = model.rts
        # self.fun = model.fun

        self.epsilonhat = self.model.get_residual()
        # print('epsilonhat',self.epsilonhat)
        # If the model is a directional distance based, set cet to CET_ADDI
        # if hasattr(self.model, 'gx'):
        #     self.model.cet = CET_ADDI
        #     self.y = np.diag(np.tensordot(
        #         self.model.y, self.model.get_gamma(), axes=([1], [1])))
        # else:

    def get_mean_of_inefficiency(self, method=RED_MOM):
        """
        Args:
            method (String, optional): RED_MOM (Method of moments) or RED_QLE (Quassi-likelihood estimation) or RED_KDE (nonparametric KDE conditional mean). Defaults to RED_MOM.
        """
        tools.assert_optimized(self.model.optimization_status)
        if method == RED_MOM:
            self.__method_of_moment(self.model.get_residual())
        elif method == RED_QLE:
            self.__quassi_likelihood(self.model.get_residual())
        elif method == RED_KDE:
            self.__gaussian_kernel_estimation(self.model.get_residual())
        else:
            raise ValueError("Undefined estimation technique.")
        return self.mu

    def get_technical_inefficiency(self, method=RED_MOM):
        """
        Args:
            method (String, optional): RED_MOM (Stevenson truncated-normal
            moment matching), RED_QLE (Stevenson truncated-normal quasi-
            likelihood), or RED_KDE (nonparametric KDE conditional mean).

        Returns firm-level conditional inefficiency u_hat = E[u_i | residual_i].
        RED_MOM and RED_QLE use the Stevenson truncated-normal JLMS formula.
        RED_KDE estimates f_u and f_v nonparametrically and computes the same
        conditional mean by numerical integration, so it no longer applies the
        Stevenson analytical JLMS formula.
        """
        tools.assert_optimized(self.model.optimization_status)
        if method == RED_KDE:
            self.get_mean_of_inefficiency(method)
            self.uhat = self.__kde_conditional_mean_inefficiency()
            return self.uhat

        if method in [RED_MOM, RED_QLE]:
            self.get_mean_of_inefficiency(method)
            sigma2 = self.sigma_u ** 2 + self.sigma_v ** 2
            sigmas = self.sigma_u * self.sigma_v / math.sqrt(sigma2)
            mu_latent = getattr(self, "mu_latent", None)
            if mu_latent is None:
                raise ValueError("Stevenson MOM/QLE requires mu_latent.")
            mus = (mu_latent * (self.sigma_v ** 2)
                   - self.residual_minus * (self.sigma_u ** 2)) / sigma2

            if hasattr(self.model, 'cet'):
                self.model.cet = self.model.cet
            else:
                self.model.cet = CET_ADDI

            z = mus / sigmas
            inverse_mills = np.exp(stats.norm.logpdf(z) - stats.norm.logcdf(z))
            jlms = sigmas * inverse_mills + mus
            self.uhat = jlms
            return self.uhat

        raise ValueError("Undefined model parameters.")

    def __jlms_efficiency_from_inefficiency(self, uhat):
        """Convert JLMS conditional inefficiency into bounded technical efficiency.

        JLMS recovers E[u | residual].  Reported technical efficiency is the
        monotone bounded transformation exp(-E[u | residual]), so it is always
        in (0, 1] for nonnegative conditional inefficiency. Tiny negative
        numerical noise is treated as zero; substantive negative values are
        reported as NaN and flagged by callers through diagnostic components.
        """
        uhat = np.asarray(uhat, dtype=float)
        tol = 1e-10
        uhat_eff = np.where((uhat < 0) & (uhat >= -tol), 0.0, uhat)
        te = np.exp(-uhat_eff)
        te = np.where(np.isfinite(uhat_eff) & (uhat_eff >= -tol), te, np.nan)
        return te

    def __standardized_distance_score(self, distance):
        """Bound a nonnegative directional distance as 1 / (1 + distance)."""
        distance = np.asarray(distance, dtype=float)
        tol = 1e-10
        distance_eff = np.where((distance < 0) & (distance >= -tol), 0.0, distance)
        score = np.divide(
            1.0,
            1.0 + distance_eff,
            out=np.full_like(distance_eff, np.nan, dtype=float),
            where=np.isfinite(distance_eff) & (distance_eff >= -tol) & np.isfinite(1.0 + distance_eff) & ((1.0 + distance_eff) > tol),
        )
        return score

    def __ddf_normalized_efficiency_from_inefficiency(self, uhat):
        """Convert DDF conditional inefficiency into bounded ratio efficiency.

        DDF ``u_hat`` is a scalar directional distance.  Its physical units
        depend on the translation anchor (basexy), NOT on any single variable.
        To form a bounded TE we normalize by the *original* (pre-translation)
        observed variable in each active direction:

        - gy active: TE_y = y / (y + u_hat * gy)
        - gx active: TE_x = (x - u_hat * gx) / x
        - gb active: TE_b = (b - u_hat * gb) / b

        Single-direction case: the corresponding TE has clear economic meaning.
        Multi-direction case: ``get_technical_efficiency()`` returns TE_binding
        = min of active component TEs; this is a binding-dimension efficiency,
        NOT an unambiguous "total TE".  Use
        ``get_technical_efficiency_components()`` for full per-direction detail.

        Important: model.x / model.y / model.b are *pre-processed* (translated)
        data that may be near zero.  We recover the original values via the
        inverse translation:  original = model_data - basexy * direction.
        """
        uhat = np.asarray(uhat, dtype=float).reshape(-1)
        tol = 1e-10
        uhat_eff = np.where((uhat < 0) & (uhat >= -tol), 0.0, uhat)
        valid_u = np.isfinite(uhat_eff) & (uhat_eff >= -tol)

        gy = np.asarray(self.gy if self.gy is not None else [], dtype=float)
        gx = np.asarray(self.gx if self.gx is not None else [], dtype=float)
        gb = np.asarray(self.gb if self.gb is not None else [], dtype=float)

        # Recover original (pre-translation) data.
        # assert_CNLSDDFweak1 applies: x_new = x + basexyb*gx, etc.
        # Inverse: original = model_data - basexyb * direction.
        base_col = np.asarray(self.basexy, dtype=float).reshape(-1, 1)

        scores = []
        score_names = []

        if gy.size and np.any(gy > 0):
            y_pre = np.asarray(self.model.y, dtype=float)
            if y_pre.ndim == 1:
                y_pre = y_pre.reshape(-1, 1)
            y_orig = y_pre + base_col * gy.reshape(1, -1)  # inverse: y_orig = y_pre - (-basexyb*gy) = y_pre + basexyb*gy... 
            # Actually assert_CNLSDDFweak1 does: y_new = y - basexyb*gy → y_orig = y_new + basexyb*gy
            # But basexy for CNLSDDFweak is basexyb; need to check sign.
            # From tools.py: y = y - basexy*gy → y_orig = y_pre + basexy*gy
            # However basexy might be basexyb (negative for gb). Let's just use abs basexy consistently.
            # Actually the translation in assert_CNLSDDFweak1 is:
            #   y_new = y - basexyb * gy
            # So: y_orig = y_new + basexyb * gy
            y_orig = y_pre + base_col * gy.reshape(1, -1)
            for idx, direction_value in enumerate(gy):
                if direction_value > 0:
                    observed = y_orig[:, idx]
                    denom = observed + uhat_eff * direction_value
                    with np.errstate(divide="ignore", invalid="ignore"):
                        score = np.divide(
                            observed, denom,
                            out=np.full_like(uhat_eff, np.nan, dtype=float),
                            where=valid_u & (observed > tol) & (denom > tol),
                        )
                    scores.append(np.clip(score, 0.0, 1.0))
                    score_names.append(f"TE_y_{idx}")

        if gx.size and np.any(gx > 0):
            x_pre = np.asarray(self.model.x, dtype=float)
            if x_pre.ndim == 1:
                x_pre = x_pre.reshape(-1, 1)
            # assert_CNLSDDFweak1: x_new = x + basexyb*gx → x_orig = x_new - basexyb*gx
            x_orig = x_pre - base_col * gx.reshape(1, -1)
            for idx, direction_value in enumerate(gx):
                if direction_value > 0:
                    observed = x_orig[:, idx]
                    projected = observed - uhat_eff * direction_value
                    with np.errstate(divide="ignore", invalid="ignore"):
                        score = np.divide(
                            projected, observed,
                            out=np.full_like(uhat_eff, np.nan, dtype=float),
                            where=valid_u & (observed > tol) & np.isfinite(projected),
                        )
                    scores.append(np.clip(score, 0.0, 1.0))
                    score_names.append(f"TE_x_{idx}")

        if gb.size and np.any(gb > 0):
            if not hasattr(self.model, "b"):
                raise ValueError("Undesirable-output direction gb requires model.b.")
            b_pre = np.asarray(self.model.b, dtype=float)
            if b_pre.ndim == 1:
                b_pre = b_pre.reshape(-1, 1)
            # assert_CNLSDDFweak1: b_new = b + basexyb*gb → b_orig = b_new - basexyb*gb
            b_orig = b_pre - base_col * gb.reshape(1, -1)
            for idx, direction_value in enumerate(gb):
                if direction_value > 0:
                    observed = b_orig[:, idx]
                    projected = observed - uhat_eff * direction_value
                    with np.errstate(divide="ignore", invalid="ignore"):
                        score = np.divide(
                            projected, observed,
                            out=np.full_like(uhat_eff, np.nan, dtype=float),
                            where=valid_u & (observed > tol) & np.isfinite(projected),
                        )
                    scores.append(np.clip(score, 0.0, 1.0))
                    score_names.append(f"TE_b_{idx}")

        if not scores:
            return self.__jlms_efficiency_from_inefficiency(uhat_eff)

        # Store per-direction scores for components()
        self._ddf_te_scores = scores
        self._ddf_te_score_names = score_names

        score_matrix = np.column_stack(scores)
        with np.errstate(all="ignore"):
            te = np.nanmin(score_matrix, axis=1)
        te = np.where(np.isfinite(te), te, np.nan)
        return te

    def __efficiency_from_inefficiency(self, uhat):
        """Dispatch to DDF normalized ratio TE or standard JLMS exp(-u_hat)."""
        if self.context.is_ddf:
            return self.__ddf_normalized_efficiency_from_inefficiency(uhat)
        return self.__jlms_efficiency_from_inefficiency(uhat)

    def __active_direction_label(self):
        """Return a compact label for active DDF direction groups."""
        gy = np.asarray(self.gy if self.gy is not None else [], dtype=float)
        gx = np.asarray(self.gx if self.gx is not None else [], dtype=float)
        gb = np.asarray(self.gb if self.gb is not None else [], dtype=float)
        active_groups = []
        if gy.size and np.any(gy > 0):
            active_groups.append("gy")
        if gx.size and np.any(gx > 0):
            active_groups.append("gx")
        if gb.size and np.any(gb > 0):
            active_groups.append("gb")
        return "+".join(active_groups) if active_groups else ""

    def get_technical_efficiency_components(self, method=RED_MOM):
        """Return JLMS efficiency and distance components for any StoNED model.

        Columns:
        - residual: model residual used by JLMS;
        - u_hat: conditional inefficiency, u_hat = E(u_i | epsilon_i);
        - directional_distance: distance-style inefficiency score, equal to
          u_hat in the normalized JLMS/DDF scale;
        - standardized_directional_distance_score: 1 / (1 + directional_distance);
        - TE: for DDF models, a real-variable-normalized direction ratio;
          otherwise JLMS technical efficiency exp(-u_hat).
        """
        tools.assert_optimized(self.model.optimization_status)

        u_hat = self.get_technical_inefficiency(method)
        directional_distance = np.asarray(u_hat, dtype=float)
        distance_score = self.__standardized_distance_score(directional_distance)
        te = self.__efficiency_from_inefficiency(u_hat)

        self.uhat = np.asarray(u_hat, dtype=float)
        self.directional_distance = directional_distance
        self.standardized_directional_distance_score = distance_score
        self.TE = te

        return pd.DataFrame({
            "direction": self.__active_direction_label(),
            "residual": np.asarray(self.epsilonhat, dtype=float),
            "u_hat": self.uhat,
            "directional_distance": directional_distance,
            "standardized_directional_distance_score": distance_score,
            "TE": te,
        })

    def get_technical_efficiency_all(self, method=RED_MOM):
        """Return JLMS efficiency and DDF distance components."""
        return self.get_technical_efficiency_components(method)

    def get_technical_efficiency(self, method=RED_MOM):
        """Return bounded TE: DDF normalized ratio or JLMS exp(-u_hat)."""
        tools.assert_optimized(self.model.optimization_status)
        self.uhat = self.get_technical_inefficiency(method)
        self.TE = self.__efficiency_from_inefficiency(self.uhat)
        return self.TE

    def __method_of_moment(self, residual):
        """Stevenson truncated-normal method of moments.

        The CNLS residual is r = v - u + E[u].  With Stevenson inefficiency
        u ~ N(mu_latent, sigma_u^2) truncated at zero, RED_MOM matches the
        second, third, and fourth cumulants of r.  This generalizes the old
        half-normal MOM; the half-normal case is the special case
        mu_latent = 0.
        """
        residual = np.asarray(residual, dtype=float)
        residual = residual[np.isfinite(residual)]
        if residual.size == 0:
            raise ValueError("Residuals must contain at least one finite value.")

        centered = residual - np.mean(residual)
        M2_mean = float(np.mean(centered ** 2))
        M3_mean = float(np.mean(centered ** 3))
        M4_mean = float(np.mean(centered ** 4))
        K4_mean = M4_mean - 3 * M2_mean ** 2

        if self.model.fun == FUN_PROD:
            sign = -1.0
            if M3_mean > max(1e-10, 1e-8 * M2_mean ** 1.5):
                raise ValueError("Production residual skewness is incompatible with nonnegative inefficiency.")
        elif self.model.fun == FUN_COST:
            sign = 1.0
            if M3_mean < -max(1e-10, 1e-8 * M2_mean ** 1.5):
                raise ValueError("Cost residual skewness is incompatible with nonnegative inefficiency.")
        else:
            raise ValueError("Undefined model parameters.")

        def __std_trunc_cumulants(eta):
            """Cumulants of Z | Z >= -eta for standard normal Z."""
            a = -eta
            q = stats.norm.cdf(eta)
            if q <= 1e-14 or not np.isfinite(q):
                return None
            phi_a = stats.norm.pdf(a)
            I0 = q
            I1 = phi_a
            I2 = a * phi_a + I0
            I3 = (a ** 2) * phi_a + 2 * I1
            I4 = (a ** 3) * phi_a + 3 * I2
            raw1 = I1 / I0
            raw2 = I2 / I0
            raw3 = I3 / I0
            raw4 = I4 / I0
            c2 = raw2 - raw1 ** 2
            c3 = raw3 - 3 * raw1 * raw2 + 2 * raw1 ** 3
            c4_central = raw4 - 4 * raw1 * raw3 + 6 * raw1 ** 2 * raw2 - 3 * raw1 ** 4
            k4 = c4_central - 3 * c2 ** 2
            if c2 <= 0 or c3 <= 0 or not np.all(np.isfinite([c2, c3, k4])):
                return None
            return c2, c3, k4

        scale2 = max(abs(M2_mean), 1e-8)
        scale3 = max(abs(M3_mean), 1e-8, 1e-6 * scale2 ** 1.5)
        scale4 = max(abs(K4_mean), 1e-8, 1e-6 * scale2 ** 2)

        half_denom = math.sqrt(2 / math.pi) * (1 - 4 / math.pi)
        if self.model.fun == FUN_PROD and M3_mean < 0:
            sigma_u0 = max((M3_mean / half_denom) ** (1 / 3), 1e-6)
        elif self.model.fun == FUN_COST and M3_mean > 0:
            sigma_u0 = max((-M3_mean / half_denom) ** (1 / 3), 1e-6)
        else:
            sigma_u0 = max(math.sqrt(M2_mean / 2), 1e-6)
        sigma_v0 = max(math.sqrt(max(M2_mean - ((math.pi - 2) / math.pi) * sigma_u0 ** 2, 1e-8)), 1e-6)

        def __moment_objective(params):
            eta, log_sigma_u, log_sigma_v = params
            sigma_u = np.exp(log_sigma_u)
            sigma_v = np.exp(log_sigma_v)
            cumulants = __std_trunc_cumulants(eta)
            if cumulants is None:
                return 1e30
            c2, c3, c4 = cumulants
            pred2 = sigma_u ** 2 * c2 + sigma_v ** 2
            pred3 = sign * sigma_u ** 3 * c3
            pred4 = sigma_u ** 4 * c4
            errors = np.array([
                (pred2 - M2_mean) / scale2,
                (pred3 - M3_mean) / scale3,
                (pred4 - K4_mean) / scale4,
            ])
            if not np.all(np.isfinite(errors)):
                return 1e30
            return float(np.sum(errors ** 2))

        starts = []
        for eta0 in [-4.0, -2.0, 0.0, 1.0, 2.0, 4.0]:
            starts.append(np.array([eta0, np.log(sigma_u0), np.log(sigma_v0)]))
            starts.append(np.array([eta0, np.log(max(math.sqrt(M2_mean), 1e-6)), np.log(max(math.sqrt(M2_mean / 2), 1e-6))]))
        bounds = [(-8, 8), (-20, 20), (-20, 20)]
        results = [
            opt.minimize(__moment_objective, start, method="L-BFGS-B", bounds=bounds)
            for start in starts
        ]
        best = min(results, key=lambda item: item.fun if np.isfinite(item.fun) else np.inf)
        if (not best.success) or (not np.isfinite(best.fun)):
            raise RuntimeError(f"Stevenson truncated-normal MOM failed: {best.message}")

        eta, log_sigma_u, log_sigma_v = best.x
        self.sigma_u = float(np.exp(log_sigma_u))
        self.sigma_v = float(np.exp(log_sigma_v))
        self.mu_latent = float(eta * self.sigma_u)
        z = self.mu_latent / self.sigma_u
        self.mu = float(self.mu_latent + self.sigma_u * np.exp(stats.norm.logpdf(z) - stats.norm.logcdf(z)))
        self.lamda = self.sigma_u / self.sigma_v
        self.moment_objective = float(best.fun)

        if self.model.fun == FUN_PROD:
            self.residual_minus = residual - self.mu
        elif self.model.fun == FUN_COST:
            self.residual_minus = -residual - self.mu

    def __quassi_likelihood(self, residual):
        residual = np.asarray(residual, dtype=float)

        if self.model.fun == FUN_PROD:
            eps_observed = residual
        elif self.model.fun == FUN_COST:
            eps_observed = -residual
        else:
            raise ValueError("Undefined model parameters.")

        eps_observed = eps_observed[np.isfinite(eps_observed)]
        if eps_observed.size == 0:
            raise ValueError("Residuals must contain at least one finite value.")

        eps_std = max(float(np.std(eps_observed, ddof=0)), 1e-8)

        def __truncated_normal_mean(mu_latent, sigma_u):
            z = mu_latent / sigma_u
            inverse_mills = np.exp(stats.norm.logpdf(z) - stats.norm.logcdf(z))
            return mu_latent + sigma_u * inverse_mills

        def __truncated_normal_loglik(params, eps):
            """Negative log likelihood for Stevenson truncated-normal SFA.

            Important residual convention:
            ``eps`` is the CNLS residual r = v - u + E[u], not the raw
            composite error e = v - u.  The likelihood below first converts
            r to the raw composite-error scale by ``centered_eps = r - E[u]``.
            If a caller already has raw e = v - u, it must use the same density
            without this centering step.

            Stevenson's inefficiency term is
            u ~ N(mu_latent, sigma_u^2) truncated at zero, where ``mu_latent``
            is the pre-truncation mean and ``self.mu`` stores E[u].
            """
            log_sigma_u, log_sigma_v, mu_latent = params
            sigma_u = np.exp(log_sigma_u)
            sigma_v = np.exp(log_sigma_v)
            sigma = np.sqrt(sigma_u ** 2 + sigma_v ** 2)
            mean_u = __truncated_normal_mean(mu_latent, sigma_u)
            # Convert CNLS residual r = v - u + E[u] to raw e = v - u.
            centered_eps = eps - mean_u

            log_denominator = stats.norm.logcdf(mu_latent / sigma_u)
            z_density = (centered_eps + mu_latent) / sigma
            z_trunc = (
                mu_latent * sigma_v ** 2 - centered_eps * sigma_u ** 2
            ) / (sigma * sigma_u * sigma_v)

            log_likelihood = (
                -np.log(sigma)
                + stats.norm.logpdf(z_density)
                + stats.norm.logcdf(z_trunc)
                - log_denominator
            )
            if not np.all(np.isfinite(log_likelihood)):
                return np.inf
            return -float(np.sum(log_likelihood))

        starts = [
            np.array([np.log(eps_std / 2), np.log(eps_std / 2), 0.0]),
            np.array([np.log(eps_std / 2), np.log(eps_std / 2), eps_std / 2]),
            np.array([np.log(eps_std / 2), np.log(eps_std / 2), -eps_std / 2]),
            np.array([np.log(eps_std), np.log(eps_std / 2), 0.0]),
            np.array([np.log(eps_std / 2), np.log(eps_std), 0.0]),
        ]
        # Put a very wide, scale-aware bound on the latent mean.  This preserves
        # Stevenson truncated-normal QLE while avoiding unbounded line-search
        # failures on nearly deterministic residual vectors.
        mu_bound = max(20.0 * eps_std, 1.0)
        bounds = [(-20, 20), (-20, 20), (-mu_bound, mu_bound)]

        results = []
        for start in starts:
            results.append(opt.minimize(
                __truncated_normal_loglik,
                start,
                args=(eps_observed,),
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 2000, "ftol": 1e-10},
            ))

        successful = [item for item in results if item.success and np.isfinite(item.fun)]
        if not successful:
            # L-BFGS-B can stop with ABNORMAL_TERMINATION_IN_LNSRCH even when
            # the likelihood is finite but numerically flat.  Powell is slower
            # but derivative-free and bounded, making it a robust fallback for
            # small residual vectors used by CQER/CER examples.
            finite_candidates = [item for item in results if np.isfinite(item.fun)]
            powell_starts = [item.x for item in finite_candidates] or starts
            for start in powell_starts:
                results.append(opt.minimize(
                    __truncated_normal_loglik,
                    start,
                    args=(eps_observed,),
                    method="Powell",
                    bounds=bounds,
                    options={"maxiter": 5000, "xtol": 1e-8, "ftol": 1e-8},
                ))
            successful = [item for item in results if item.success and np.isfinite(item.fun)]

        finite_results = [item for item in results if np.isfinite(item.fun)]
        candidates = successful if successful else finite_results
        if not candidates:
            message = results[0].message if results else "no optimizer result"
            raise RuntimeError(f"Stevenson truncated-normal QLE failed: {message}")

        best = min(candidates, key=lambda item: item.fun)

        log_sigma_u, log_sigma_v, mu_latent = best.x
        self.sigma_u = float(np.exp(log_sigma_u))
        self.sigma_v = float(np.exp(log_sigma_v))
        self.mu_latent = float(mu_latent)
        self.mu = float(__truncated_normal_mean(self.mu_latent, self.sigma_u))
        self.lamda = self.sigma_u / self.sigma_v
        self.log_likelihood = -float(best.fun)

        if self.model.fun == FUN_PROD:
            self.residual_minus = residual - self.mu
        elif self.model.fun == FUN_COST:
            self.residual_minus = -residual - self.mu

    def __kde_oriented_residual(self, residual):
        """Orient CNLS residuals to eps_obs = v - u + E[u]."""
        residual = np.asarray(residual, dtype=float)
        if self.model.fun == FUN_PROD:
            return residual
        if self.model.fun == FUN_COST:
            return -residual
        raise ValueError("Undefined model parameters.")

    def __safe_gaussian_kde(self, sample):
        """Build a Gaussian KDE; add deterministic jitter only for degenerate samples."""
        sample = np.asarray(sample, dtype=float).reshape(-1)
        sample = sample[np.isfinite(sample)]
        if sample.size < 2:
            return None
        std = float(np.std(sample, ddof=1)) if sample.size > 1 else 0.0
        if std <= 1e-12:
            spread = max(abs(float(sample[0])), 1.0) * 1e-6
            sample = sample + np.linspace(-spread, spread, sample.size)
        return stats.gaussian_kde(sample)

    def __kde_reflected_u_pdf(self, kde_u, u_grid):
        """Boundary-corrected KDE for u >= 0 using reflection at zero."""
        if kde_u is None:
            return np.zeros_like(u_grid, dtype=float)
        density = kde_u(u_grid) + kde_u(-u_grid)
        return np.where(np.isfinite(density) & (density > 0), density, 0.0)

    def __kde_initial_inefficiency_sample(self, eps_observed):
        """Construct a nonnegative one-sided inefficiency seed for KDE.

        The seed estimates f_u and f_v nonparametrically.  The returned
        firm-level u_hat is not this seed; it is the numerical conditional
        mean E[u|eps].
        """
        eps_observed = np.asarray(eps_observed, dtype=float).reshape(-1)
        finite = eps_observed[np.isfinite(eps_observed)]
        if finite.size < 3:
            raise ValueError("RED_KDE requires at least three finite residuals.")
        eps_std = max(float(np.std(finite, ddof=1)), 1e-8)
        frontier = float(np.quantile(finite, 0.95))
        u_seed = np.maximum(frontier - finite, 0.0)
        if np.std(u_seed, ddof=1) <= 1e-12 or np.mean(u_seed) <= 1e-12:
            frontier = float(np.max(finite))
            u_seed = np.maximum(frontier - finite, 0.0)
        if np.std(u_seed, ddof=1) <= 1e-12:
            ranks = stats.rankdata(-finite, method="average") / finite.size
            u_seed = eps_std * np.maximum(ranks - np.min(ranks), 0.0)
        return u_seed

    def __kde_conditional_mean_inefficiency(self):
        """Compute RED_KDE u_hat = E[u | eps] by numerical integration.

        RED_KDE does not impose Stevenson truncated-normal inefficiency. Given
        KDE estimates of f_u(u) on u >= 0 and f_v(v), it evaluates
            E[u|eps] = int u f_v(eps+u) f_u(u) du / int f_v(eps+u) f_u(u) du.
        """
        eps_raw = np.asarray(self.residual_minus, dtype=float).reshape(-1)
        u_sample = np.asarray(self.kde_u_sample, dtype=float).reshape(-1)
        v_sample = np.asarray(self.kde_v_sample, dtype=float).reshape(-1)
        finite_u = u_sample[np.isfinite(u_sample) & (u_sample >= 0)]
        finite_v = v_sample[np.isfinite(v_sample)]
        if finite_u.size < 3 or finite_v.size < 3:
            raise ValueError("RED_KDE conditional mean requires finite u and v KDE samples.")
        if np.nanmax(finite_u) <= 1e-12:
            return np.zeros_like(eps_raw, dtype=float)

        kde_u = self.__safe_gaussian_kde(finite_u)
        kde_v = self.__safe_gaussian_kde(finite_v)
        if kde_u is None or kde_v is None:
            raise ValueError("RED_KDE could not construct KDE densities for u and v.")

        u_upper = max(
            float(np.quantile(finite_u, 0.995)) * 1.5,
            float(np.max(finite_u)) * 1.1,
            float(np.std(finite_u, ddof=1)) * 8.0,
            1e-8,
        )
        grid_count = int(min(1200, max(300, finite_u.size * 4)))
        u_grid = np.linspace(0.0, u_upper, grid_count)
        fu = self.__kde_reflected_u_pdf(kde_u, u_grid)
        if not np.any(fu > 0):
            raise ValueError("RED_KDE estimated a zero inefficiency density on the integration grid.")

        out = np.full_like(eps_raw, np.nan, dtype=float)
        for idx, eps_value in enumerate(eps_raw):
            if not np.isfinite(eps_value):
                continue
            fv = kde_v(eps_value + u_grid)
            integrand = fu * fv
            integrand = np.where(np.isfinite(integrand) & (integrand > 0), integrand, 0.0)
            den = np.trapz(integrand, u_grid)
            if not np.isfinite(den) or den <= 1e-300:
                continue
            num = np.trapz(u_grid * integrand, u_grid)
            out[idx] = max(float(num / den), 0.0)

        self.kde_integration_grid = u_grid
        return out

    def __gaussian_kernel_estimation(self, residual):
        """Nonparametric KDE residual decomposition for RED_KDE.

        This path does not estimate Stevenson truncated-normal parameters. It
        orients the residual to eps_obs = v - u + E[u], builds a nonnegative
        one-sided seed for u, estimates f_u and f_v by KDE, and lets
        get_technical_inefficiency(RED_KDE) recover firm-level u_hat by
        numerical conditional expectation.
        """
        eps_observed = self.__kde_oriented_residual(residual)
        eps_observed = np.asarray(eps_observed, dtype=float).reshape(-1)
        finite_mask = np.isfinite(eps_observed)
        finite_eps = eps_observed[finite_mask]
        if finite_eps.size < 3:
            raise ValueError("KDE residual decomposition requires at least three finite residuals.")

        u_seed = self.__kde_initial_inefficiency_sample(finite_eps)
        mu = float(np.mean(u_seed))
        eps_raw_finite = finite_eps - mu
        v_seed = eps_raw_finite + u_seed
        v_seed = v_seed - np.mean(v_seed)

        residual_minus = np.full_like(eps_observed, np.nan, dtype=float)
        residual_minus[finite_mask] = eps_observed[finite_mask] - mu

        self.mu = mu
        self.mu_latent = np.nan
        self.sigma_u = float(np.std(u_seed, ddof=1)) if u_seed.size > 1 else 0.0
        self.sigma_v = float(np.std(v_seed, ddof=1)) if v_seed.size > 1 else 0.0
        self.lamda = self.sigma_u / self.sigma_v if self.sigma_v > 1e-12 else np.inf
        self.kde_u_sample = np.maximum(u_seed, 0.0)
        self.kde_v_sample = v_seed
        self.kde_method = "nonparametric_conditional_mean"
        self.residual_minus = residual_minus


    # def get_SDFDDFhat(self, method=RED_KDE):
    #     """
    #     Args:
    #         method (String, optional): RED_MOM (Method of moments) or RED_QLE (Quassi-likelihood estimation). Defaults to RED_MOM.
    #
    #     Calculate the StoNED frontier
    #     """
    #     tools.assert_optimized(self.model.optimization_status)
    #     self.get_mean_of_inefficiency(method)
    #
    #     model2 = CNLSSDFDDF_FRONTIER.CNLSDDF(self.data, sent=self.sent,
    #                    muhat=self.mu, epsilonhat=self.epsilonhat, z=self.z, gy=self.gy, gx=self.gx, fun = self.fun, rts=self.rts)
    #     model2.optimize(solver="mosek")
    #     SDFDDF = model2.get_obj_minus_ddf()
    #     return  SDFDDF


    def richardson_lucy_blind_corrected(self, method):
        # Legacy Richardson-Lucy deconvolution helper. Current RED_KDE uses
        # nonparametric KDE conditional-mean integration and does not call this.
        """
        修正后的Richardson-Lucy Blind Deconvolution算法，用于从残差epsilon中估计u和v。

        参数：
        -----------
        epsilon : 1D numpy数组
            残差，用于估计效率u和噪声v。
        kernel_size : int
            噪声核v的大小。
        max_m : int
            最大的盲迭代次数。
        max_j : int
            每个盲迭代中的RL迭代次数。
        tol : float
            收敛容差。
        verbose : bool
            如果为True，则打印收敛信息。

        返回：
        --------
        u_final : 1D numpy数组
            估计的企业特定效率。
        v_final : 1D numpy数组
            估计的噪声核。
        """
        kernel_size = 5
        max_m = 500
        max_j = 500
        tol = 1e-16
        verbose = True
        # 确保epsilon是1D numpy数组
        tools.assert_optimized(self.model.optimization_status)
        self.get_mean_of_inefficiency(method)
        # print("????????????????")
        epsilon = self.model.get_residual() + self.mu
        epsilon = np.asarray(epsilon, dtype=np.float64).flatten()
        N = len(epsilon)

        # 选择一个足够大的M以确保epsilon_shifted为非负
        M = max(0, -np.min(epsilon)) + 1.0
        epsilon_shifted = epsilon + M

        # 初始化u和v
        u = epsilon_shifted.copy()  # 按描述初始化u为调整后的残差
        v = np.ones(kernel_size, dtype=np.float64) / kernel_size  # 初始化v为全1向量并归一化

        # 小常数防止除以零
        eps = 1e-12

        for m in range(max_m):
            if verbose:
                print(f"盲迭代 {m+1}/{max_m}")

            # 备份当前u和v以检查收敛
            u_prev = u.copy()
            v_prev = v.copy()

            # Step 2.1: 更新v
            for j in range(max_j):
                # 计算当前的预测epsilon
                u_convolved_v = convolve(u, v, mode='same') + eps

                # 计算E
                E = epsilon_shifted / u_convolved_v

                # 互相关u和E来更新v
                correlation = correlate(u, E, mode='full')

                # 截取与kernel_size对齐的部分
                center = len(correlation) // 2
                start = center - kernel_size // 2
                end = start + kernel_size
                correlation = correlation[start:end]

                # 更新v
                v *= correlation

                # 防止v中出现全部为零的情况
                sum_v = np.sum(v)
                if sum_v > 0:
                    v /= sum_v  # 归一化v
                else:
                    v = np.ones(kernel_size, dtype=np.float64) / kernel_size  # 重置v并归一化

                if verbose and j == 0:
                    print(f"  更新v的第 {j+1}/{max_j} 次RL迭代")

            # Step 2.2: 更新u
            for j in range(max_j):
                # 重新计算预测epsilon
                u_convolved_v = convolve(u, v, mode='same') + eps

                # 计算E
                E = epsilon_shifted / u_convolved_v

                # 卷积v和E来更新u
                convolution = convolve(E, v[::-1], mode='same')  # 旋转v以实现互相关

                # 更新u
                u *= convolution

                # 确保u非负
                u = np.maximum(u, 0)

                if verbose and j == 0:
                    print(f"  更新u的第 {j+1}/{max_j} 次RL迭代")

            # 检查收敛
            delta_u = np.linalg.norm(u - u_prev) / (np.linalg.norm(u_prev) + eps)
            delta_v = np.linalg.norm(v - v_prev) / (np.linalg.norm(v_prev) + eps)
            if verbose:
                print(f"  收敛信息: Δu = {delta_u:.6e}, Δv = {delta_v:.6e}")
            if delta_u < tol and delta_v < tol:
                if verbose:
                    print("  已达到收敛条件，提前终止迭代。")
                break

        # 最终的u需要减去M
        u_final = u - M

        # 确保u_final非负
        u_final = np.maximum(u_final, 0)

        return u_final
