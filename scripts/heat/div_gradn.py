#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
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

import numpy as np
from soma import aims

heatmap_volume = aims.read("heat.nii.gz")
header = heatmap_volume.header()
voxel_size_x, voxel_size_y, voxel_size_z = header["voxel_size"][:3]

heat = np.asarray(heatmap_volume)

class HalfcentredGradients:
    def __init__(self, array):
        self.ppp = array[1:, 1:, 1:]
        self.ppm = array[1:, 1:, :-1]
        self.pmp = array[1:, :-1, 1:]
        self.pmm = array[1:, :-1, :-1]
        self.mpp = array[:-1, 1:, 1:]
        self.mpm = array[:-1, 1:, :-1]
        self.mmp = array[:-1, :-1, 1:]
        self.mmm = array[:-1, :-1, :-1]

    def gradx(self):
        return (0.25 / voxel_size_x) * (self.ppp + self.pmp + self.ppm + self.pmm
                                        - (self.mpp + self.mmp + self.mpm + self.mmm))
    def grady(self):
        return (0.25 / voxel_size_y) * (self.ppp + self.mpp + self.ppm + self.mpm
                                        - (self.pmp + self.mmp + self.pmm + self.mmm))
    def gradz(self):
        return (0.25 / voxel_size_z) * (self.ppp + self.mpp + self.pmp + self.mmp
                                        - (self.ppm + self.mpm + self.pmm + self.mmm))
    def gradxyz(self):
        return self.gradx(), self.grady(), self.gradz()

grad = HalfcentredGradients(heat)
gradx, grady, gradz = grad.gradxyz()
gradn = np.sqrt(gradx ** 2 + grady ** 2 + gradz ** 2)
with np.errstate(divide="ignore", invalid="ignore"):
    gradx /= gradn
    grady /= gradn
    gradz /= gradn

grad = HalfcentredGradients(gradx)
ggradx = grad.gradx()
grad = HalfcentredGradients(grady)
ggrady = grad.grady()
grad = HalfcentredGradients(gradz)
ggradz = grad.gradz()

div = ggradx + ggrady + ggradz

div_volume = aims.Volume(heatmap_volume)
div_volume.fill(float("NaN"))
div_array = np.asarray(div_volume)

div_array[1:-1, 1:-1, 1:-1] = div

# This is needed since AIMS does not silently replace NaN values in input files
# by zeros anymore. Without it, ylAdvectTubes ends up interpolating NaN values,
# which ends badly.
# TODO: remove this horrible hack!!!
div_array[np.isnan(div_array)] = 0

aims.write(div_volume, "heat_div_gradn.nii.gz")
