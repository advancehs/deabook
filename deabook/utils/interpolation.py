import numpy as np


def interpolation(alpha, beta, x_test, fun=None):
    x_arr = np.asarray(x_test, dtype=float)
    beta_arr = np.asarray(beta, dtype=float)
    alpha_arr = np.asarray(alpha, dtype=float).reshape(-1)
    if x_arr.ndim == 1:
        x_arr = x_arr.reshape(1, -1)
    if beta_arr.ndim == 1:
        beta_arr = beta_arr.reshape(1, -1)
    # Evaluate each test point against all local hyperplanes and take the
    # concave envelope convention used by production CNLS.
    vals = x_arr @ beta_arr.T + alpha_arr
    return np.min(vals, axis=1)
