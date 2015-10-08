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

import os.path

import numpy
from soma import aims

from highres_cortex.cortex_topo import CSF_LABEL, CORTEX_LABEL, WHITE_LABEL

def difference_in_cortex(image1, image2, classif):
    return numpy.abs(image2 - image1)[classif == CORTEX_LABEL]

def max_cortex_difference_files(result_file, reference_file, classif):
    result_vol = aims.read(result_file)
    result = numpy.asarray(result_vol)

    reference_vol = aims.read(reference_file)
    reference = numpy.asarray(reference_vol)

    difference = difference_in_cortex(result, reference, classif)

    return numpy.max(difference)


def compare_in_directory(dir):
    classif_vol = aims.read(os.path.join(dir, "classif.nii.gz"))
    classif = numpy.asarray(classif_vol)

    largest_voxel_size = max(classif_vol.getVoxelSize()[:3])
    real_thickness = aims.read(os.path.join(dir, "reference_thickness.nii.gz")).value(0, 0, 0)
    print("For reference: 1 voxel = {0:.2f}mm = {1:.1f}% of cortical thickness"
          .format(largest_voxel_size,
                  100 * largest_voxel_size / real_thickness))

    print("Euclidean max difference (advection): {0:.1f}%"
          .format(100 * max_cortex_difference_files(
                os.path.join(dir, "laplace-euclidean", "pial-fraction.nii.gz"),
                os.path.join(dir, "reference_euclidean.nii.gz"),
                classif)))
    print("Euclidean max difference (upwinding): {0:.1f}%"
          .format(100 * max_cortex_difference_files(
                os.path.join(dir, "upwind-euclidean", "pial-fraction.nii.gz"),
                os.path.join(dir, "reference_euclidean.nii.gz"),
                classif)))

    abs_thickness_diff = max_cortex_difference_files(
                os.path.join(dir, "laplace-euclidean", "total-length.nii.gz"),
                os.path.join(dir, "reference_thickness.nii.gz"),
                classif)
    print("Thickness max difference (advection): {0:.2f}mm ({1:.1f}%)"
          .format(abs_thickness_diff,
                  100 * abs_thickness_diff / real_thickness))
    abs_thickness_diff = max_cortex_difference_files(
                os.path.join(dir, "upwind-euclidean", "total-length.nii.gz"),
                os.path.join(dir, "reference_thickness.nii.gz"),
                classif)
    print("Thickness max difference (upwinding): {0:.2f}mm ({1:.1f}%)"
          .format(abs_thickness_diff,
                  100 * abs_thickness_diff / real_thickness))

    print("Laplacian value max difference: {0:.1f}%"
          .format(100 * max_cortex_difference_files(
                os.path.join(dir, "heat", "heat.nii.gz"),
                os.path.join(dir, "reference_laplacian.nii.gz"),
                classif)))
    print("Equivolumic max difference (advection): {0:.1f}%"
          .format(100 * max_cortex_difference_files(
                os.path.join(dir, "isovolume", "pial-volume-fraction.nii.gz"),
                os.path.join(dir, "reference_equivolumic.nii.gz"),
                classif)))

if __name__ == "__main__":
    compare_in_directory(".")
