
import numpy as np
import pandas as pd
from pyomo.environ import Constraint, Set


def cg_pairs(*matrices):
    """Return (i, h) pairs selected by one or more active matrices."""
    if not matrices:
        return []
    shape = None
    selected = None
    for matrix in matrices:
        if matrix is None:
            continue
        arr = np.asarray(matrix)
        if arr.size == 0:
            continue
        if shape is None:
            shape = arr.shape
            selected = np.zeros(shape, dtype=bool)
        selected |= arr.astype(float) != 0
    if selected is None:
        return []
    return [(int(i), int(h)) for i, h in zip(*np.where(selected)) if int(i) != int(h)]


def replace_constraint(obj, name, pairs, rule, doc=None):
    """Replace a full pairwise Pyomo constraint with a sparse pair set."""
    model = obj.__model__
    if hasattr(model, name):
        model.del_component(getattr(model, name))
    set_name = '_' + name + '_pairs'
    if hasattr(model, set_name):
        model.del_component(getattr(model, set_name))
    setattr(model, set_name, Set(dimen=2, initialize=list(pairs)))
    pair_set = getattr(model, set_name)
    setattr(model, name, Constraint(pair_set, rule=rule, doc=doc or name))


def make_wp_dataframe(y, x, b, z=None):
    y_arr = np.asarray(y, dtype=float).reshape(-1)
    x_arr = np.asarray(x, dtype=float)
    b_arr = np.asarray(b, dtype=float)
    if x_arr.ndim == 1:
        x_arr = x_arr.reshape(-1, 1)
    if b_arr.ndim == 1:
        b_arr = b_arr.reshape(-1, 1)
    data = pd.DataFrame({f'x{j}': x_arr[:, j] for j in range(x_arr.shape[1])})
    data['y0'] = y_arr
    for j in range(b_arr.shape[1]):
        data[f'b{j}'] = b_arr[:, j]
    z_name = None
    if z is not None:
        z_arr = np.asarray(z, dtype=float)
        if z_arr.ndim == 1:
            z_arr = z_arr.reshape(-1, 1)
        z_cols = []
        for j in range(z_arr.shape[1]):
            col = f'z{j}'
            data[col] = z_arr[:, j]
            z_cols.append(col)
        z_name = ' '.join(z_cols)
    sent = ' '.join(f'x{j}' for j in range(x_arr.shape[1])) + '=y0:' + ' '.join(f'b{j}' for j in range(b_arr.shape[1]))
    return data, sent, z_name
