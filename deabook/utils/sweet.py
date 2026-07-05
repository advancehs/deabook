# import dependencies
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay


def _as_frame_for_sweet(x, prefix="i"):
    if isinstance(x, (pd.DataFrame, pd.Series)):
        return x
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return pd.DataFrame(arr, index=list(range(arr.shape[0])))


def sweet(x, xref=None):
    x = _as_frame_for_sweet(x)
    xref = x if xref is None else _as_frame_for_sweet(xref)
    n, m = x.shape[0], xref.shape[0]
    cutactive = np.zeros((n, m))

    # Always include diagonal/common-index self constraints when available.
    common = [idx for idx in x.index if idx in set(xref.index)]
    for idx in common:
        try:
            cutactive[list(x.index).index(idx), list(xref.index).index(idx)] = 1
        except ValueError:
            pass

    arr = np.asarray(xref, dtype=float)
    if m <= 1:
        return pd.DataFrame(cutactive, index=x.index, columns=xref.index)

    try:
        # Delaunay can fail for low-rank or tiny samples; nearest-neighbor fallback below is then used.
        if arr.shape[0] > arr.shape[1] + 1:
            tri = Delaunay(arr)
            indptr, indices = tri.vertex_neighbor_vertices
            for j in range(m):
                neighbors = indices[indptr[j]:indptr[j + 1]]
                for h in neighbors:
                    if j < n and h < m:
                        cutactive[j, h] = 1
    except Exception:
        pass

    # Fallback: add each evaluated point's nearest reference point if no off-diagonal was found.
    if np.sum(cutactive) <= len(common):
        x_arr = np.asarray(x, dtype=float)
        for i in range(n):
            d = np.sum((arr - x_arr[i, :]) ** 2, axis=1)
            h = int(np.argmin(d))
            cutactive[i, h] = 1

    return pd.DataFrame(cutactive, index=x.index, columns=xref.index)
