# import dependencies
import numpy as np
import pandas
import pandas as pd
import math
import scipy.stats as stats
import scipy.optimize as opt
from .utils import tools
from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RED_MOM, RED_QLE, RED_KDE
from pyomo.environ import   exp
from math import sqrt, pi, log
from scipy.signal import convolve, correlate

class StoNED:
    """Stochastic nonparametric envelopment of data (StoNED)
    """

    def __init__(self, model):
        """StoNED
        model: The input model for residual decomposition
        """

        # print(model.__class__.__name__) # CNLSSD


        self.model = model
        if model.__class__.__name__ in ["CNLSSDweak","CNLSDDFweak",]:
            self.basexy = model.basexyb
            self.gb = model.gb

        elif model.__class__.__name__ in ["CNLSSD","CNLSDDF"]:
            self.basexy = model.basexy

        elif model.__class__.__name__ in ["CNLSSDweakmeta", "CNLSDDFweakmeta"]:
            self.basexy = model.basexyb
            # print('aaabasexyb', self.basexy)

            self.gb = model.gb
            self.gddf_er = model.gddf_er
            self.gddf = model.gddf

            self.gresidual =model.gresidual
            self.basexy_old = model.basexyb_old

        self.y = model.y
        # self.data = model.data
        # self.sent = model.sent
        # self.z = model.z
        self.gy = model.gy
        self.gx = model.gx
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
            method (String, optional): RED_MOM (Method of moments) or RED_QLE (Quassi-likelihood estimation) or RED_KDE (Kernel deconvolution estimation). Defaults to RED_MOM.
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
            method (String, optional): RED_MOM (Method of moments) or RED_QLE (Quassi-likelihood estimation). Defaults to RED_MOM.

        calculate sigma_u, sigma_v, mu, and epsilon value
        """
        tools.assert_optimized(self.model.optimization_status)
        if method == RED_QLE:
            self.get_mean_of_inefficiency(method)
            # print("mean_ie",self.get_mean_of_inefficiency(method))
            sigmas = self.sigma_u * self.sigma_v / math.sqrt(self.sigma_u ** 2 + self.sigma_v ** 2)
            mus = (self.mu * (self.sigma_v**2)-self.residual_minus * (self.sigma_u**2)) / (self.sigma_u ** 2 + self.sigma_v ** 2)

            if hasattr(self.model, 'cet'):
                self.model.cet = self.model.cet
            else:
                self.model.cet = CET_ADDI
            if self.model.fun == FUN_PROD:
                if self.model.cet == CET_ADDI:

                    jlms = sigmas * (stats.norm.pdf(mus / sigmas)) / (stats.norm.cdf(mus / sigmas)) + mus
                    # bc = np.exp(-mus + 0.5 * (sigmas) ** 2) * (stats.norm.cdf((mus / sigmas) - sigmas) / stats.norm.cdf(mus / sigmas))

                    return jlms
                elif self.model.cet == CET_MULT:

                    # bc = np.exp(-mus + 0.5 * (sigmas) ** 2) * (
                    #         stats.norm.cdf((mus / sigmas) - sigmas) / stats.norm.cdf(mus / sigmas))  ## 这个计算的是技术效率！！！！！

                    bc = sigmas * (stats.norm.pdf(mus / sigmas)) / (stats.norm.cdf(mus / sigmas)) + mus

                    return bc
        elif method == RED_KDE:
            self.get_mean_of_inefficiency(method)
            u = self.richardson_lucy_blind_corrected(method)
            return u

        raise ValueError("Undefined model parameters.")

    def get_technical_efficiency_ratio(self, method=RED_MOM):  ## for DDF
        """
        Args:
            method (String, optional): RED_MOM (Method of moments) or RED_QLE (Quassi-likelihood estimation). Defaults to RED_MOM.

        calculate sigma_u, sigma_v, mu, and epsilon value
        """
        tools.assert_optimized(self.model.optimization_status)
        self.ddfhat = self.get_technical_inefficiency(method)

        if self.model.__class__.__name__ in ["CNLSDDFweak","CNLSDDF"]:

            if sum(self.gy) >= 1:
                self.TE = self.basexy/(self.basexy+self.ddfhat) ## self.ddfhat >0
                return self.TE

            elif sum(self.gx) >= 1:

                # XXX = (-1 * np.asarray(self.basexy)  + self.epsilonhat)
                # for i, row in enumerate(XXX):
                #     if row < 0:
                #         print(f"XXX[{i}] = {row} 小于0，将其替换为011", -1 * np.asarray(self.basexy)[i],  self.epsilonhat[i],ddfhat[i])
                
                self.TE = (-1*np.asarray(self.basexy) + self.epsilonhat) /(-1*np.asarray(self.basexy) + self.epsilonhat + self.ddfhat)

                # self.TE = (-1*np.asarray(self.basexy) -self.ddfhat) /(-1*np.asarray(self.basexy))
                return self.TE

            elif sum(self.gb) >= 1:
                # self.TE = (-1*np.asarray(self.basexy) -self.ddfhat) /(-1*np.asarray(self.basexy))
                self.TE = (-1*np.asarray(self.basexy) + self.epsilonhat) /(-1*np.asarray(self.basexy) + self.epsilonhat + self.ddfhat)

                return self.TE


        elif self.model.__class__.__name__ in ["CNLSSDweakmeta","CNLSDDFweakmeta"]:

            if sum(self.gy) >= 1:
                self.TGR = self.basexy/(self.basexy+self.ddfhat)
                df = pd.DataFrame({
                    'TGR': self.TGR,
                    # 'TGDR': [1-i for i in self.TGR],
                    'GTE': self.gddf_er,
                    'MTE': [ (i*j) for i, j in zip(self.TGR, self.gddf_er)],
                })

                return df.reset_index(drop=True)

            elif sum(self.gx) >= 1:

                # XXX = (-1 * np.asarray(self.basexy) - self.gddf + self.epsilonhat)
                # for i, row in enumerate(XXX):
                #     if row < 0:
                #             print(f"XXX[{i}] = {row} 小于0，将其替换为0", -1 * np.asarray(self.basexy)[i], self.gddf[i], self.epsilonhat[i],ddfhat[i])

                self.TGR = (-1 * np.asarray(self.basexy) - self.gddf + self.epsilonhat
                             )/ (-1 * np.asarray(self.basexy) - self.gddf + self.epsilonhat + ddfhat )
                # self.TGR = (-1*np.asarray(self.basexy) -self.ddfhat) /(-1*np.asarray(self.basexy))

                df = pd.DataFrame({
                    'TGR': self.TGR,
                    'GTE': self.gddf_er,
                    'MTE': [ (i*j) for i, j in zip(self.TGR, self.gddf_er)],
                })
                return df.reset_index(drop=True)

            elif sum(self.gb) >= 1:
                self.TGR = (-1 * np.asarray(self.basexy) -  self.gddf- self.epsilonhat
                             )/ (-1 * np.asarray(self.basexy) -  self.gddf- self.epsilonhat + ddfhat )
                # self.TGR = (-1*np.asarray(self.basexy) -self.ddfhat) /(-1*np.asarray(self.basexy))

                df = pd.DataFrame({
                    'TGR': self.TGR,
                    'GTE': self.gddf_er,
                    'MTE': [ (i*j) for i, j in zip(self.TGR, self.gddf_er)],
                })
                return df.reset_index(drop=True)

        else:
            raise ValueError("Undefined model parameters.")




    def get_technical_efficiency(self, method=RED_MOM):
        """
        Args:
            method (String, optional): RED_MOM (Method of moments) or RED_QLE (Quassi-likelihood estimation). Defaults to RED_MOM.

        calculate sigma_u, sigma_v, mu, and epsilon value
        """
        tools.assert_optimized(self.model.optimization_status)
        if self.model.cet == CET_ADDI:

            te_jlms = np.exp(-self.get_technical_inefficiency(method))
            return te_jlms

        elif self.model.cet == CET_MULT:

            # te_bc = self.get_technical_inefficiency(method)
            te_bc =  np.exp(-self.get_technical_inefficiency(method))

            return te_bc

        raise ValueError("Undefined model parameters.")

    def __method_of_moment(self, residual):
        """Method of moment"""
        M2 = (residual - np.mean(residual)) ** 2
        M3 = (residual - np.mean(residual)) ** 3

        M2_mean = np.mean(M2, axis=0)
        M3_mean = np.mean(M3, axis=0)

        if self.model.fun == FUN_PROD:
            if M3_mean > 0:
                M3_mean = 0.0
            self.sigma_u = (M3_mean / ((2 / math.pi) ** (1 / 2) *
                                       (1 - 4 / math.pi))) ** (1 / 3)

        elif self.model.fun == FUN_COST:
            if M3_mean < 0:
                M3_mean = 0.00001
            self.sigma_u = (-M3_mean / ((2 / math.pi) ** (1 / 2) *
                                        (1 - 4 / math.pi))) ** (1 / 3)

        else:
            raise ValueError("Undefined model parameters.")

        self.sigma_v = (M2_mean - ((math.pi - 2) / math.pi) * self.sigma_u ** 2) ** (1 / 2)
        self.mu = (self.sigma_u ** 2 * 2 / math.pi) ** (1 / 2)

        # print("mu",self.mu)
        if self.model.fun == FUN_PROD:
            self.residual_minus = residual - self.mu
        else:
            self.residual_minus = residual + self.mu


    def __quassi_likelihood(self, residual):
        def __quassi_likelihood_estimation(lamda, eps):
            """ This function computes the negative of the log likelihood function
            given parameter (lambda) and residual (eps).

            Args:
                lamda (float): signal-to-noise ratio
                eps (list): values of the residual

            Returns:
                float: -logl, negative value of log likelihood
            """
            # sigma Eq. (3.26) in Johnson and Kuosmanen (2015)
            sigma = np.sqrt(
                np.mean(eps ** 2) / (1 - 2 * lamda ** 2 / (pi * (1 + lamda ** 2))))

            # bias adjusted residuals Eq. (3.25)
            # mean
            mu = sqrt(2 / pi) * sigma * lamda / sqrt(1 + lamda ** 2)

            # adj. res.
            epsilon = eps - mu

            # log-likelihood function Eq. (3.24)
            pn = stats.norm.cdf(-epsilon * lamda / sigma)
            return -(-len(epsilon) * log(sigma) + np.sum(np.log(pn)) -
                     0.5 * np.sum(epsilon ** 2) / sigma ** 2)

        if self.model.fun == FUN_PROD:
            lamda = opt.minimize(__quassi_likelihood_estimation,
                                 1.0,
                                 residual,
                                 method='BFGS').x[0]
        elif self.model.fun == FUN_COST:
            lamda = opt.minimize(__quassi_likelihood_estimation,
                                 1.0,
                                 -residual,
                                 method='BFGS').x[0]
        else:
            # TODO(error/warning handling): Raise error while undefined fun
            return False
        # use estimate of lambda to calculate sigma Eq. (3.26) in Johnson and Kuosmanen (2015)
        sigma = sqrt(np.mean(residual ** 2) /
                     (1 - (2 * lamda ** 2) / (pi * (1 + lamda ** 2))))

        # calculate bias correction
        # (unconditional) mean
        self.mu = sqrt(2) * sigma * lamda / sqrt(pi * (1 + lamda ** 2))

        # calculate sigma.u and sigma.v
        self.sigma_v = (sigma ** 2 / (1 + lamda ** 2)) ** (1 / 2)
        self.sigma_u = self.sigma_v * lamda

        if self.model.fun == FUN_PROD:
            self.residual_minus = residual - self.mu
        elif self.model.fun == FUN_COST:
            self.residual_minus = residual + self.mu

    def __gaussian_kernel_estimation(self, residual):
        def __gaussian_kernel_estimator(g):
            """Gaussian kernel estimator"""
            return (1 / math.sqrt(2 * math.pi)) * np.exp(-0.5 * g ** 2)

        residual=-residual
        x = np.array(residual)

        # choose a bandwidth (rule-of-thumb, Eq. (3.29) in Silverman (1986))
        if np.std(x, ddof=1) < stats.iqr(x, interpolation='midpoint'):
            estimated_sigma = np.std(x, ddof=1)
        else:
            estimated_sigma = stats.iqr(x, interpolation='midpoint')
        h = 1.06 * estimated_sigma * len(self.y) ** (-1 / 5)

        # kernel matrix
        kernel_matrix = np.zeros((len(self.y), len(self.y)))
        for i in range(len(self.y)):
            kernel_matrix[i] = np.array([
                __gaussian_kernel_estimator(g=((x[i] - x[j]) / h)) /
                (len(self.y) * h) for j in range(len(self.y))
            ])

        # kernel density value
        kernel_density_value = np.sum(kernel_matrix, axis=0)

        # unconditional expected inefficiency mu
        derivative = np.zeros(len(self.y))
        for i in range(len(self.y) - 1):
            derivative[i +
                       1] = 0.2 * (kernel_density_value[i + 1] -
                                   kernel_density_value[i]) / (x[i + 1] - x[i])


        # expected inefficiency mu
        self.mu = np.max(derivative)
        # print("selfmu",self.mu)
        if self.model.fun == FUN_PROD:
            self.residual_minus = residual + self.mu
            # print(1111111111111,self.residual_minus)

        elif self.model.fun == FUN_COST:
            self.residual_minus = residual - self.mu




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
