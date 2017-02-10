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

import math
import os.path
import sys

import numpy
from soma import aims

from highres_cortex.cortex_topo import CSF_LABEL, CORTEX_LABEL, WHITE_LABEL

def difference_from_files(result_file, reference_file, classif_array):
    result_vol = aims.read(result_file)
    result = numpy.asarray(result_vol)

    reference_vol = aims.read(reference_file)
    reference = numpy.asarray(reference_vol)

    difference = numpy.ma.masked_array(result - reference,
                                       mask=(classif_array != CORTEX_LABEL))

    return difference

def symmetrize_value_range(value_range, center):
    breadth = max(value_range[1] - center, center - value_range[0])
    return (center - breadth, center + breadth)

def scatter_plot_files(result_file, reference_file, classif,
                       value_range=None, range_centre=None, ax=None):
    result_vol = aims.read(result_file)
    result_values = numpy.asarray(result_vol)[classif == CORTEX_LABEL]

    reference_vol = aims.read(reference_file)
    reference_values = numpy.asarray(reference_vol)[classif == CORTEX_LABEL]

    if ax is None:
        import matplotlib.pyplot
        ax = matplotlib.pyplot.figure().add_subplot(111)
    ax.scatter(reference_values, result_values)
    if value_range is None:
        value_range = (min(reference_values.min(), result_values.min()),
                       max(reference_values.max(), result_values.max()))
    if range_centre is not None:
        value_range = symmetrize_value_range(value_range, range_centre)
    ax.plot(value_range, value_range)
    ax.set_xlim(value_range)
    ax.set_ylim(value_range)
    ax.set_aspect("equal")


class ResultComparator:
    reference_file = {
        os.path.join("heat", "heat.nii.gz"):
            "reference_laplacian.nii.gz",
        os.path.join("heat", "heat_div_gradn.nii.gz"):
            "reference_curvature.nii.gz",
        os.path.join("dist", "distCSF.nii.gz"):
            "reference_distCSF.nii.gz",
        os.path.join("dist", "distwhite.nii.gz"):
            "reference_distwhite.nii.gz",
        os.path.join("dist", "dist_sum.nii.gz"):
            "reference_thickness.nii.gz",
        os.path.join("laplace-euclidean", "total-length.nii.gz"):
            "reference_thickness.nii.gz",
        os.path.join("laplace-euclidean", "pial-fraction.nii.gz"):
            "reference_euclidean.nii.gz",
        os.path.join("isovolume", "pial-volume-fraction.nii.gz"):
            "reference_equivolumic.nii.gz",
        os.path.join("upwind-euclidean", "total-length.nii.gz"):
            "reference_thickness.nii.gz",
        os.path.join("upwind-euclidean", "pial-fraction.nii.gz"):
            "reference_euclidean.nii.gz",
        os.path.join("upwind-equivolume", "corrected-pial-volume-fraction.nii.gz"):
            "reference_equivolumic.nii.gz",
        os.path.join("CBS", "Equivolumic", "_surf_thickness.nii.gz"):
            "reference_thickness.nii.gz",
        os.path.join("CBS", "Equidistant", "inverted_layering.nii.gz"):
            "reference_euclidean.nii.gz",
        os.path.join("CBS", "Equivolumic", "inverted_layering.nii.gz"):
            "reference_equivolumic.nii.gz",
    }

    dimension = {
        os.path.join("heat", "heat.nii.gz"): "%",
        os.path.join("heat", "heat_div_gradn.nii.gz"): "mm^{-1}",
        os.path.join("dist", "distwhite.nii.gz"): "mm",
        os.path.join("dist", "distCSF.nii.gz"): "mm",
        os.path.join("dist", "dist_sum.nii.gz"): "mm",
        os.path.join("laplace-euclidean", "total-length.nii.gz"): "mm",
        os.path.join("laplace-euclidean", "pial-fraction.nii.gz"): "%",
        os.path.join("isovolume", "pial-volume-fraction.nii.gz"): "%",
        os.path.join("upwind-euclidean", "total-length.nii.gz"): "mm",
        os.path.join("upwind-euclidean", "pial-fraction.nii.gz"): "%",
        os.path.join("upwind-equivolume", "corrected-pial-volume-fraction.nii.gz"): "%",
        os.path.join("CBS", "Equivolumic", "_surf_thickness.nii.gz"): "mm",
        os.path.join("CBS", "Equidistant", "inverted_layering.nii.gz"): "%",
        os.path.join("CBS", "Equivolumic", "inverted_layering.nii.gz"): "%",
    }

    def __init__(self, dir):
        self._dir = dir
        path = self._make_subpath
        self._classif_vol = aims.read(path("classif.nii.gz"))
        self._classif = numpy.asarray(self._classif_vol)
        self._voxel_size = tuple(self._classif_vol.getVoxelSize()[:3])
        thickness_vol = aims.read(path("reference_thickness.nii.gz"))
        self._thickness = thickness_vol.value(0, 0, 0)

    def _make_subpath(self, *components):
        return os.path.join(self._dir, *components)

    def scatter_plot_file(self, result_file, ax=None, *args, **kwargs):
        path = self._make_subpath

        if ax is None:
            import matplotlib.pyplot
            ax = matplotlib.pyplot.figure().add_subplot(111)

        reference_file = self.reference_file[result_file]
        scatter_plot_files(path(result_file), path(reference_file),
                           self._classif, ax=ax, *args, **kwargs)
        ax.set_title(result_file)

    def plot_compare_all(self, fig=None):
        path = self._make_subpath

        if fig is None:
            import matplotlib.pyplot
            fig = matplotlib.pyplot.figure()

        if os.path.isdir(path("CBS")):
            include_CBS = True
            num_lines = 3
        else:
            include_CBS = False
            num_lines = 2

        #ax = fig.add_subplot(num_lines, 4, 1)
        #self.scatter_plot_file(os.path.join("heat", "heat.nii.gz"),
        #                       value_range=(0, 1), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 1)
        self.scatter_plot_file(os.path.join("dist", "distCSF.nii.gz"),
                               value_range=(0, 1.1*self._thickness), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 5)
        self.scatter_plot_file(os.path.join("dist", "dist_sum.nii.gz"),
                               value_range=(0, 2 * self._thickness), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 2)
        self.scatter_plot_file(os.path.join("laplace-euclidean", "total-length.nii.gz"),
                               value_range=(0, 2 * self._thickness), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 6)
        self.scatter_plot_file(os.path.join("upwind-euclidean", "total-length.nii.gz"),
                               value_range=(0, 2 * self._thickness), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 3)
        self.scatter_plot_file(os.path.join("laplace-euclidean", "pial-fraction.nii.gz"),
                               value_range=(0, 1), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 7)
        self.scatter_plot_file(os.path.join("upwind-euclidean", "pial-fraction.nii.gz"),
                               value_range=(0, 1), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 4)
        self.scatter_plot_file(os.path.join("isovolume", "pial-volume-fraction.nii.gz"),
                               value_range=(0, 1), ax=ax)

        ax = fig.add_subplot(num_lines, 4, 8)
        try:
            self.scatter_plot_file(os.path.join("upwind-equivolume", "corrected-pial-volume-fraction.nii.gz"),
                                value_range=(0, 1), ax=ax)
        except IOError:
            pass

        if include_CBS:
            try:
                ax = fig.add_subplot(num_lines, 4, 10)
                self.scatter_plot_file(os.path.join("CBS", "Equivolumic", "_surf_thickness.nii.gz"),
                                       value_range=(0, 2 * self._thickness), ax=ax)
            except IOError:
                pass
            try:
                ax = fig.add_subplot(num_lines, 4, 11)
                self.scatter_plot_file(os.path.join("CBS", "Equidistant", "inverted_layering.nii.gz"),
                                       value_range=(0, 1), ax=ax)
            except IOError:
                pass
            try:
                ax = fig.add_subplot(num_lines, 4, 12)
                self.scatter_plot_file(os.path.join("CBS", "Equivolumic", "inverted_layering.nii.gz"),
                                       value_range=(0, 1), ax=ax)
            except IOError:
                pass


    def compare_files(self, result_file, reference_file=None):
        path = self._make_subpath

        if reference_file is None:
            reference_file = self.reference_file[os.path.normpath(result_file)]
        diff = difference_from_files(
            path(result_file), path(reference_file),
            self._classif)

        rms_error = math.sqrt((diff ** 2).mean())
        bias = diff.mean()

        return (rms_error, bias)

    @classmethod
    def comparison_to_text(self, rms_error, bias, result_file=None):
        try:
            dimension = self.dimension[result_file]
        except KeyError:
            dimension = ""
        if dimension == "%":
            return ("RMS error = {0:.1f}%, bias = {1:.1f}%"
                    .format(100 * rms_error, 100 * bias))
        elif dimension:
            return ("RMS error = {0:.3f} {unit}, bias = {1:.3f} {unit}"
                    .format(rms_error, bias, unit=dimension))
        else:
            return ("RMS error = {0:.3g}, bias = {1:.3g}"
                    .format(rms_error, bias))

    def text_compare_files(self, result_file, reference_file=None):
        rms_error, bias = self.compare_files(result_file, reference_file)
        return self.comparison_to_text(rms_error, bias, result_file)

    def ensure_max_rms_error(self, result_file, max_rms_error,
                             reference_file=None):
        rms_error, bias = self.compare_files(result_file, reference_file)
        text = "{0}: {1}".format(
            result_file,
            self.comparison_to_text(rms_error, bias, result_file))

        if rms_error <= max_rms_error:
            text += " (RMS error <= {0})".format(max_rms_error)
        else:
            text += " (RMS error > {0}) <== ERROR".format(max_rms_error)

        print(text)
        return rms_error <= max_rms_error

    def ensure_max_rms_errors(self, list_of_tests):
        success = True
        for test_params in list_of_tests:
            ret = self.ensure_max_rms_error(*test_params)
            if not ret:
                success = False
        return success

    def text_compare_all(self):
        print("voxel size = {0:.2f} mm = {1:.1f}% of cortical thickness"
              .format(max(self._voxel_size),
                      100 * max(self._voxel_size) / self._thickness))

        for result_file in sorted(self.reference_file.iterkeys()):
            sys.stdout.write("{0}: ".format(result_file))
            try:
                print(self.text_compare_files(result_file))
            except IOError:
                print("-- cannot read file")
                continue


if __name__ == "__main__":
    comparator = ResultComparator(".")
    comparator.text_compare_all()
    comparator.plot_compare_all()
    import matplotlib.pyplot
    matplotlib.pyplot.show()
