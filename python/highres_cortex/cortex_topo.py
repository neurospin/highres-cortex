# Copyright CEA (2014).
# Copyright Université Paris XI (2014).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# This file is part of highres-cortex, a collection of software designed
# to process high-resolution magnetic resonance images of the cerebral
# cortex.
#
# This software is governed by the CeCILL licence under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# licence as circulated by CEA, CNRS and INRIA at the following URL:
# <http://www.cecill.info/>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the licence, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of scientific
# software, that may mean that it is complicated to manipulate, and that
# also therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL licence and that you accept its terms.

"""Tools for manipulating the topology of the cortex.
"""

import numpy as np
from soma import aims
from soma import aimsalgo


FLT_MAX = 3.4028234663852886e+38
NaN = float("NaN")

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
    fm = aims.FastMarching("6")
    fm.setVerbose(verbose)
    dist = fm.doit(classif_vol, propagation_labels, seed_labels)
    np_dist = np.asarray(dist)
    mask = (np_dist == FLT_MAX)
    #np_dist[mask] = NaN
    np_classif[np_dist == 0] = border_label

    # The connectivity does not seem to matter here
    fm = aims.FastMarching("6")
    fm.setVerbose(verbose)
    dist_neg = fm.doit(classif_vol, seed_labels,
                       list(propagation_labels) + [border_label])
    np_dist_neg = np.asarray(dist_neg)
    #np_dist_neg[np_dist_neg == FLT_MAX] = NaN
    np_dist[mask] = -np_dist_neg[mask]

    return dist
