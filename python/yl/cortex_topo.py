"""Tools for manipulating the topology of the cortex.
"""

import numpy as np
from soma import aims
from soma import aimsalgo


FLT_MAX = 3.4028234663852886e+38


def fastmarching_negative(classif_vol,
                          propagation_labels=[1],
                          seed_labels=[0],
                          border_label=5,
                          verbose=True):
    """Distance to the seed as well as negative distance within the seed.

    Caution: the input volume classif_vol is modified: the value border_label
    is given to voxels at the border (seed side).
    """
    # Positive distances are propagated from the seed (seed_labels). The
    # distance is considered to be zero at the seed's side of the boundary:
    # therefore, this layer needs to be added to the cortex to make the seed
    # for negative distances. 6-connectivity is used for the boundary so that
    # the distance gradient is continuous across this boundary.
    np_classif = np.asarray(classif_vol)
    fm = aims.FastMarching_Volume_S16_6c()
    fm.setVerbose(verbose)
    dist = fm.doit(classif_vol, propagation_labels, seed_labels)
    np_dist = np.asarray(dist)
    np_classif[np_dist == 0] = border_label

    # The connectivity does not seem to matter here
    fm = aims.FastMarching_Volume_S16_6c()
    fm.setVerbose(verbose)
    dist_neg = fm.doit(classif_vol, seed_labels,
                       list(propagation_labels) + [border_label])
    np_dist_neg = np.asarray(dist_neg)
    mask = (np_dist == FLT_MAX)
    np_dist[mask] = -np_dist_neg[mask]

    if len(seed_labels) == 1:
        print(np.any(np_classif[np_dist_neg == 0] == seed_labels[0]))

    return dist
