# import dependencies
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from .tools import trans_list, to_2d_list


def sweet(x):
    """Sweet spot approach

    Args:
        x (DataFrame): input  pandas.

    Returns:
        list: active concavity constraint.
    """

    # transform data
    df = np.asmatrix(x)

    # calculate distance matrix
    distance = cdist(df, df)
    distance[np.diag_indices_from(distance)] = np.nan

    # calculate distance cut
    distcut = np.asmatrix(np.nanpercentile(distance, 3, axis=0))
    print(distcut)
    # find concavity constraint in sweet spot
    distance = np.where(np.isnan(distance), 0, distance)

    cutactive = np.zeros((distance.shape[0], distance.shape[1]))
    for i in range(distance.shape[0]):
        for j in range(distance.shape[1]):
            # print("1",distcut[:, i],distcut[:, j])
            if distance[i, j] <= distcut[:, i]:
                cutactive[i, j] = 1
    cutactive2 = pd.DataFrame(cutactive,index=x.index,columns=x.index)

    return cutactive2


def sweetref(x,xref):
    """Sweet spot approach

    Args:
        x (DataFrame): input  pandas.
        xref (DataFrame): input reference pandas.

    Returns:
        list: active concavity constraint.
    """
    # transform data
    df = np.asmatrix(x)
    df2 = np.asmatrix(xref)

    # calculate distance matrix
    distance = cdist(df, df2)
    # distance[np.diag_indices_from(distance)] = np.nan

    # calculate distance cut
    distcut = np.asmatrix(np.nanpercentile(distance, 3, axis=0))

    # find concavity constraint in sweet spot
    distance = np.where(np.isnan(distance), 0, distance)

    cutactive = np.zeros((distance.shape[0], distance.shape[1]))
    for i in range(distance.shape[0]):
        for j in range(distance.shape[1]):
            if distance[i, j] <= distcut[:, i]:
                cutactive[i, j] = 1
    cutactive2 = pd.DataFrame(cutactive,index=x.index,columns=xref.index)
    return cutactive2
