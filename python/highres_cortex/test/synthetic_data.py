# -*- coding: utf-8 -*-
# Copyright Télécom ParisTech (2015).
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

import math

from soma import aims
import numpy

def _convert_to_float_triple(param):
    try:
        value = float(param)
        return (value, value, value)
    except TypeError:
        pass

    assert len(param) == 3
    return tuple(float(item) for item in param)

def make_cortex_sphere_classif(inner_radius, outer_radius,
                               voxel_size, margin=None):
    inner_radius = _convert_to_float_triple(inner_radius)
    outer_radius = _convert_to_float_triple(outer_radius)
    voxel_size = _convert_to_float_triple(voxel_size)
    if margin is None:
        margin = max(voxel_size)
    margin = _convert_to_float_triple(margin)

    for i in range(3):
        assert outer_radius[i] > inner_radius[i] > 0

    size = [int(math.ceil(2 * (ax_radius + ax_margin) / ax_voxel_size))
            for ax_radius, ax_margin, ax_voxel_size
            in zip(outer_radius, margin, voxel_size)]

    classif_volume = aims.Volume(size[0], size[1], size[2], dtype="S16")
    classif_volume.header()["voxel_size"] = [
        voxel_size[0], voxel_size[1], voxel_size[2], 1]
    emplace_cortex_sphere_classif(classif_volume,
                                  inner_radius, outer_radius, margin)
    return classif_volume

def make_centred_coord_grids(classif_volume):
    size = (classif_volume.getSizeX(),
            classif_volume.getSizeY(),
            classif_volume.getSizeZ())
    voxel_size = classif_volume.header()["voxel_size"][:3]

    s = [numpy.linspace(-ax_voxel_size * (ax_size // 2),
                        ax_voxel_size * (ax_size // 2),
                        ax_size)
         for ax_voxel_size, ax_size in zip(voxel_size, size)]
    return numpy.ix_(*s)


def emplace_cortex_sphere_classif(classif_volume,
                                  inner_radius, outer_radius, margin):
    inner_radius = _convert_to_float_triple(inner_radius)
    outer_radius = _convert_to_float_triple(outer_radius)
    margin = _convert_to_float_triple(margin)

    size = (classif_volume.getSizeX(),
            classif_volume.getSizeY(),
            classif_volume.getSizeZ())
    voxel_size = classif_volume.header()["voxel_size"][:3]
    for i in range(3):
        assert outer_radius[i] > inner_radius[i] > 0
        assert size[i] * voxel_size[i] >= 2 * outer_radius[i]

    grids = make_centred_coord_grids(classif_volume)

    np_classif = numpy.asarray(classif_volume)
    classif_volume.fill(0)
    np_classif[sum(grid ** 2 / radius ** 2
                   for radius, grid in zip(outer_radius, grids)) < 1] = 100
    np_classif[sum(grid ** 2 / radius ** 2
                   for radius, grid in zip(inner_radius, grids)) < 1] = 200
