# -*- coding: utf-8 -*-
# Copyright Forschungszentrum Jülich GmbH (2017).
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
                               voxel_size, margin=None,
                               noise=None, sigma=None):
    voxel_size = _convert_to_float_triple(voxel_size)
    if margin is None:
        margin = max(voxel_size)
    margin = _convert_to_float_triple(margin)

    assert outer_radius > inner_radius > 0

    size = [int(math.ceil(2 * (outer_radius + ax_margin) / ax_voxel_size))
            for ax_margin, ax_voxel_size
            in zip(margin, voxel_size)]

    classif_aimsdata = aims.AimsData(size[0], size[1], size[2], dtype="S16")
    classif_aimsdata.setSizeXYZT(*voxel_size)
    classif_volume = classif_aimsdata.volume()
    emplace_cortex_sphere_classif(classif_volume,
                                  inner_radius, outer_radius,
                                  margin=margin, noise=noise, sigma=sigma)
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
                                  inner_radius, outer_radius,
                                  margin, noise, sigma):
    margin = _convert_to_float_triple(margin)

    size = (classif_volume.getSizeX(),
            classif_volume.getSizeY(),
            classif_volume.getSizeZ())
    voxel_size = classif_volume.header()["voxel_size"][:3]
    assert outer_radius > inner_radius > 0
    for i in range(3):
        assert size[i] * voxel_size[i] >= 2 * outer_radius

    grids = make_centred_coord_grids(classif_volume)

    np_classif = numpy.asarray(classif_volume)
    classif_volume.fill(0)

    distance_to_centre = sum(grid ** 2 for grid in grids)
    if noise:
        distance_to_centre += make_noise_array(noise, sigma, size)
    np_classif[distance_to_centre < outer_radius ** 2] = 100
    np_classif[distance_to_centre < inner_radius ** 2] = 200


def make_noise_array(noise, sigma, size):
    import scipy.ndimage.filters
    arr = noise * numpy.random.randn(*size)
    scipy.ndimage.filters.gaussian_filter(arr, sigma, output=arr)
    return arr


def _make_similar_volume(data_array, ref):
    volume = aims.Volume(data_array)
    volume.header().update(ref.header())
    return volume


def make_sphere_and_reference_result(inner_radius, outer_radius, voxel_size,
                                     margin=None, noise=None, sigma=None):
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    assert outer_radius > inner_radius > 0
    voxel_size = _convert_to_float_triple(voxel_size)

    thickness = outer_radius - inner_radius

    classif_volume = make_cortex_sphere_classif(inner_radius, outer_radius,
                                                voxel_size, margin=margin,
                                                noise=noise, sigma=sigma)
    grids = make_centred_coord_grids(classif_volume)
    distance_to_centre = numpy.sqrt(sum(grid ** 2 for grid in grids))

    distance_to_white = distance_to_centre - inner_radius
    distance_to_CSF = outer_radius - distance_to_centre

    euclidean_metric = numpy.clip((outer_radius - distance_to_centre)
                                  / (outer_radius - inner_radius),
                                  0, 1)

    with numpy.errstate(divide="ignore"):
        laplacian_value = numpy.clip(
            inner_radius / (outer_radius - inner_radius) *
            (outer_radius / distance_to_centre - 1),
            0, 1)
        curvature = - 2 / distance_to_centre

    equivolumic_metric = numpy.clip(
        (outer_radius ** 3 - distance_to_centre ** 3) /
        (outer_radius ** 3 - inner_radius ** 3),
        0, 1)

    return (classif_volume,
            _make_similar_volume(distance_to_white, ref=classif_volume),
            _make_similar_volume(distance_to_CSF, ref=classif_volume),
            _make_similar_volume(euclidean_metric, ref=classif_volume),
            _make_similar_volume(laplacian_value, ref=classif_volume),
            _make_similar_volume(curvature, ref=classif_volume),
            _make_similar_volume(equivolumic_metric, ref=classif_volume))


def write_sphere_and_reference_result(inner_radius, outer_radius, voxel_size,
                                      dir=".", margin=None,
                                      noise=None, sigma=None):
    if not os.path.isdir(dir):
        os.makedirs(dir)
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    assert outer_radius > inner_radius > 0
    voxel_size = _convert_to_float_triple(voxel_size)

    (classif,
     distance_to_white, distance_to_CSF,
     euclidean_metric,
     laplacian_value,
     curvature,
     equivolumic_metric) = (
         make_sphere_and_reference_result(
             inner_radius, outer_radius, voxel_size, margin=margin,
             noise=noise, sigma=sigma))

    np_thickness = (outer_radius - inner_radius) * (
        numpy.ones_like(distance_to_white))
    thickness = _make_similar_volume(np_thickness, ref=distance_to_white)

    aims.write(classif,
               os.path.join(dir, "classif.nii.gz"))
    aims.write(distance_to_white,
               os.path.join(dir, "reference_distwhite.nii.gz"))
    aims.write(distance_to_CSF,
               os.path.join(dir, "reference_distCSF.nii.gz"))
    aims.write(thickness,
               os.path.join(dir, "reference_thickness.nii.gz"))
    aims.write(euclidean_metric,
               os.path.join(dir, "reference_euclidean.nii.gz"))
    aims.write(laplacian_value,
               os.path.join(dir, "reference_laplacian.nii.gz"))
    aims.write(curvature,
               os.path.join(dir, "reference_curvature.nii.gz"))
    aims.write(equivolumic_metric,
               os.path.join(dir, "reference_equivolumic.nii.gz"))


if __name__ == "__main__":
    import os
    import shutil
    import argparse
    parser = argparse.ArgumentParser(
        description="Write a synthetic sphere and reference results "
                    "to the current directory.")
    parser.add_argument("inner_radius", type=float)
    parser.add_argument("thickness", type=float)
    parser.add_argument("voxel_size", type=float)
    parser.add_argument("--margin", type=float, default=None)
    parser.add_argument("--noise", type=float, default=None)
    parser.add_argument("--sigma", type=float, default=None)
    args = parser.parse_args()
    write_sphere_and_reference_result(
        args.inner_radius,
        args.inner_radius + args.thickness,
        args.voxel_size,
        margin=args.margin,
        noise=args.noise,
        sigma=args.sigma)
