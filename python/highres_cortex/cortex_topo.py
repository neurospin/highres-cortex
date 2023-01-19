# -*- coding: utf-8 -*-
# Copyright Télécom ParisTech (2015).
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

from __future__ import absolute_import, division, print_function

import os
import tempfile
import shutil
import soma.subprocess as subprocess

import numpy as np
from soma import aims
from soma import aimsalgo


FLT_MAX = 3.4028234663852886e+38
NaN = np.nan

CSF_LABEL = 0
CORTEX_LABEL = 100
WHITE_LABEL = 200


def signed_distance(classif_vol,
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
    # np_dist[mask] = NaN
    np_classif[np_dist == 0] = border_label

    # The connectivity does not seem to matter here
    fm = aims.FastMarching("6")
    fm.setVerbose(verbose)
    dist_neg = fm.doit(classif_vol, seed_labels,
                       list(propagation_labels) + [border_label])
    np_dist_neg = np.asarray(dist_neg)
    # np_dist_neg[np_dist_neg == FLT_MAX] = NaN
    np_dist[mask] = -np_dist_neg[mask]

    return dist


def fix_cortex_topology(input_classif, filling_size=2., fclosing=10.):
    """Fix the topology of a cortical segmentation.

    The topology of a hollow sphere is imposed onto a voxelwise segmentation of
    the cortex, which consists of the following labels:

    Label 0 (`CSF_LABEL`)
        Outside of the cortex, corresponding to the cerebrospinal fluid,
        defined in 26-connectivity.

    Label 100 (`CORTEX_LABEL`)
        The cortex itself, defined using 6-connectivity.

    Label 200 (`WHITE_LABEL`)
        Inside of the cortex, corresponds to the white matter,
        defined in 26-connectivity.


    Parameters
    ----------

    classif: aims.Volume
        The input voxelwise classification.

    filling_size: float
        The size, in millimetres, of the largest holes in either cortical
        boundary that will be filled. This must be kept smaller than the
        smallest cortical thickness in the image (see `Method` below for a more
        precise description of this parameter). The default value is 2 mm,
        which is appropriate for a human brain.

    fclosing: float
        The radius, in millimetres, of the morphological closing which is used
        by VipHomotopic in Cortical surface mode to retrieve the brain's outer
        envelope. The default value, 10 mm, is appropriate for a human brain.

    Returns
    -------

    The topology-corrected voxelwise classification is returned in an
    `aims.Volume_S16`.

    Raises
    ------

    OSError
        This function throws ``OSError`` if ``VipHomotopic`` cannot be found
        or executed.

    soma.subprocess.CalledProcessError
        This exception can occur if ``VipHomotopic``, which is in charge of the
        homotopic morphological operations, terminates with an error.


    Environment
    -----------

    This function needs the ``VipHomotopic`` command from the Morphologist
    image segmentation pipeline to reside in the ``PATH``. Note that the
    original ``VipHomotopic`` has hard-coded limits on the number of iterations
    for morphological operations, which may be exceeded when working on
    high-resolution (sub-millimetre) images.

    Input/output
    ------------

    The command ``VipHomotopic``, which is used to perform the homotopic
    morphological operations, reports progress on stdout/stderr.

    Images are passed to ``VipHomotopic`` using files under a temporary
    directory allocated with `tempfile.mkdtemp`.

    Method
    ------

    The topology correction is done in two main steps:

    1. A topologically spherical bounding box of the brain is computed and
    dilated towards the inside until it reaches the white matter. This
    retrieves a topologically correct object which fits the grey--white
    boundary.

    2. The previous object is eroded from the inside in the region where it
    overlaps with the cortex. This retrieves a topologically correct pial
    boundary.

    Each of these main steps is performed in two sub-steps: first the homotopic
    morpholological operation is performed until a boundary which is dilated by
    `filling_size`, then to the original boundary. This guides the front
    propagation, in order to prevent the formation of spurious strands.

    This method will change voxels from the cortex class to either the white
    matter or the CSF class, as needed to ensure the topology.

    Note that the output is not a deterministic function of the input, because
    the homotopic operations use a pseudo-random order for the front
    propagation.
    """
    fclosing = float(fclosing)
    assert fclosing >= 0
    filling_size = float(filling_size)
    assert filling_size >= 0

    # VipHomotopic only works with 16-bit signed integer voxels.
    conv = aims.ShallowConverter(intype=input_classif, outtype="Volume_S16")
    classif = conv(input_classif)

    tmp_classif = _prepare_classif_for_VipHomotopic_Cortical(classif,
                                                             filling_size)

    tmp_dir = None
    try:
        tmp_dir = tempfile.mkdtemp(prefix="highres-cortex.")
        aims.write(tmp_classif, os.path.join(tmp_dir, "tmp_classif.nii.gz"))
        del tmp_classif
        with open(os.path.join(tmp_dir, "fake.han"), "w") as f:
            f.write("sequence: unknown\n"
                    "gray: mean: 120 sigma: 10\n"
                    "white: mean: 433 sigma: 10\n")
        # VipHomotopic in Cortical surface mode retrieves a spherical
        # grey--white boundary by iteratively eroding the bounding box of the
        # cortex in a homotopic manner. It will proceed in two steps, first
        # stopping at STEP1_FRONT_BARRIER, and finally at WHITE_LABEL.
        subprocess.check_call(["VipHomotopic", "-mode", "C",
                               "-input", "tmp_classif.nii.gz",
                               "-classif", "tmp_classif.nii.gz",
                               "-hana", "fake.han",
                               "-fclosing", repr(fclosing),
                               "-output", "cortex.nii.gz"], cwd=tmp_dir)

        aims.write(classif, os.path.join(tmp_dir, "classif.nii.gz"))

        # First boundary to guide VipHomotopic (prevent leaking through holes
        # in sulci).
        aimsdata_classif = aims.Volume_S16(classif.getSize, [1, 1, 1])
        aimsdata_classif[:] = classif[:]
        # Restore the header (in particular the voxel_size), which may not have
        # been copied in the constructor because a border is requested.
        aimsdata_classif.header().update(classif.header())
        eroded = aimsalgo.AimsMorphoErosion(aimsdata_classif, filling_size)
        del classif, aimsdata_classif
        aims.write(eroded, os.path.join(tmp_dir, "eroded.nii.gz"))
        del eroded

        # The spherical grey--white boundary is dilated in a homotopic manner
        # until the border of eroded_classif is reached.
        subprocess.check_call(["VipHomotopic", "-mode", "H",
                               "-input", "eroded.nii.gz",
                               "-cortex", "cortex.nii.gz",
                               "-fclosing", "0",
                               "-output", "bigsulci.nii.gz"], cwd=tmp_dir)
        subprocess.check_call(["VipHomotopic", "-mode", "H",
                               "-input", "classif.nii.gz",
                               "-cortex", "bigsulci.nii.gz",
                               "-fclosing", "0",
                               "-output", "pial_surface.nii.gz"], cwd=tmp_dir)

        cortex = aims.read(os.path.join(tmp_dir, "cortex.nii.gz"))
        pial_surface = aims.read(os.path.join(tmp_dir, "pial_surface.nii.gz"))
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    array_cortex = np.asarray(cortex)
    array_pial_surface = np.asarray(pial_surface)

    array_cortex[array_cortex == 0] = 200
    array_cortex[array_cortex == 255] = 100
    array_cortex[array_pial_surface != 0] = 0
    return cortex


def _prepare_classif_for_VipHomotopic_Cortical(classif, filling_size):
    array_classif = np.asarray(classif)

    # Dilate the cortex by 1 voxel, otherwise the topologically spherical front
    # initialized by VipHomotopic (Cortical surface mode) is at the inside of
    # the border, which means that the final segmentation will be one voxel
    # off. The default -fclosing parameter prevents this from happening in
    # sulci by closing them, but it still happens at the top of gyri. Setting a
    # fake voxel size of 1 mm is a trick for expressing the dilation in terms
    # of voxels.
    #
    # The 1-voxel border is necessary for AimsMorpho{Dilation,Erosion}.
    aimsdata_classif = aims.Volume_S16(classif.getSize(), [1, 1, 1])
    aimsdata_classif[:] = classif[:]
    saved_voxel_size = classif.header()["voxel_size"][:3]
    aimsdata_classif.setVoxelSize(1, 1, 1)
    dilated = aimsalgo.AimsMorphoDilation(aimsdata_classif, 1)
    # Restore the voxel size in case the header is shared with the aims.Volume
    # that aimsdata_classif was created from (classif). BUG: restoring the
    # value like this is not thread-safe!
    aimsdata_classif.setVoxelSize(saved_voxel_size)
    del aimsdata_classif
    array_dilated = dilated.np

    tmp_classif = aims.Volume(classif)
    array_tmp_classif = np.asarray(tmp_classif)
    array_tmp_classif[np.logical_and(array_tmp_classif == CSF_LABEL,
                                     array_dilated != 0)] = CORTEX_LABEL
    del dilated, array_dilated

    # This almost-white region will serve as a first boundary in VipHomotopic,
    # before the final dilation towards the white matter. This helps restore
    # the continuity of the white matter so that the homotopic criterion
    # chooses the right connections (i.e. do not create spurious strands or
    # planes through the cortex).
    white = aims.Volume_S16(classif.gteSize(), [1, 1, 1])
    white[:] = classif[:]
    # Restore the header (in particular the voxel_size), which may not have
    # been copied in the constructor because a border is requested.
    white.header().update(classif.header())
    array_white = white.np
    array_white[array_white != WHITE_LABEL] = CSF_LABEL
    dilated_white = aimsalgo.AimsMorphoDilation(white, filling_size)
    del white, array_white
    array_dilated_white = dilated_white.np
    STEP1_FRONT_BARRIER = 199
    array_tmp_classif[array_dilated_white != 0] = STEP1_FRONT_BARRIER
    del dilated_white, array_dilated_white
    array_tmp_classif[array_classif == WHITE_LABEL] = WHITE_LABEL

    return tmp_classif
