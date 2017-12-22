# -*- coding: utf-8 -*-
# Copyright Forschungszentrum Jülich GmbH (2017).
#
# Author: Yann Leprince <yann.leprince@ylep.fr>.
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


"""Elementary processes to be used within the CAPSUL pipelining system."""


import capsul.api
import capsul.process.xml
from traits.api import Bool, Enum, File, Float, Int, Str, Undefined


VOLUME_EXTENSIONS = ['.nii.gz', '.vimg', '.vinfo', '.vhdr', '.img', '.hdr',
                     '.v', '.i', '.mnc', '.mnc.gz', '.nii', '.jpg', '.gif',
                     '.png', '.mng', '.bmp', '.pbm', '.pgm', '.ppm', '.xbm',
                     '.xpm', '.tiff', '.tif', '.ima', '.dim', '']


class Laplacian(capsul.api.Process):
    """Solve the Laplacian model in the cortex"""

    classif = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")
    precision = Float(
        0.001, output=False, optional=True,
        desc="target maximum relative error in first-order finite differences")
    typical_cortical_thickness = Float(
        3, output=False, optional=True,
        desc="typical thickness of the cortex (mm), used for accelerating "
        "convergence")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    laplace_field = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output pseudo-temperature field (from 0 in CSF to 1 in the "
        "white matter)"
    )

    def get_commandline(self):
        return [
            "ylLaplacian",
            "--classif", self.classif,
            "--output", self.laplace_field,
            "--precision", repr(self.precision),
            "--typical-cortical-thickness",
            repr(self.typical_cortical_thickness),
            "--verbose", str(self.verbosity)]


class IsoCurvature(capsul.api.Process):
    """Compute the curvature of isosurfaces"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image volume (scalar field)")
    # modes mean, geom, pri1, pri2 are currently unimplemented
    mode = Enum(
        "sum", values=("sum",), output=False, optional=True,
        desc="type of curvature to compute")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output image volume containing the curvature of the isosurfaces "
        "of the input field"
    )

    def get_commandline(self):
        return [
            "ylIsoCurvature",
            "--input", self.input_image,
            "--mode", self.mode,
            "--output", self.output_image,
            "--verbose", str(self.verbosity)]


class RemoveNaN(capsul.api.Process):
    """Remove NaN values from an image"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    value = Float(
        0, output=False, optional=True,
        desc="replacement value")
    percentage = Bool(
        True, output=False, optional=True,
        desc="interpret value as a percentage of the image intensity range")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output image"
    )

    def get_commandline(self):
        return [
            "AimsRemoveNaN",
            "-i", self.input_image,
            "-np", str(self.percentage),
            "--value", repr(self.value),
            "-o", self.output_image,
            "--verbose", str(self.verbosity)]


class MedianFilter(capsul.api.Process):
    """Median filter smoothing"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    x_size = Int(3, output=False, optional=True, desc="X size of the filter mask")
    y_size = Int(3, output=False, optional=True, desc="Y size of the filter mask")
    z_size = Int(3, output=False, optional=True, desc="Z size of the filter mask")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="median-filtered image"
    )

    def get_commandline(self):
        return [
            "AimsMedianSmoothing",
            "--verbose", "0",
            "--input", self.input_image,
            "--dx", str(self.x_size),
            "--dy", str(self.y_size),
            "--dz", str(self.z_size),
            "--output", self.output_image]


class BinarizeCortex(capsul.api.Process):
    """Extract a binary image (0/1) of the cortex"""

    classif = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="binary image of the cortex (1 in the cortex, 0 elsewhere)"
    )

    def get_commandline(self):
        return [
            "AimsThreshold",
            "--verbose", "0",
            "-b",
            "--fg", "1",
            "-m", "eq",
            "-t", "100",
            "--input", self.classif,
            "--output", self.output_image]


class AdvectTubesAlongGradient(capsul.api.Process):
    """Advect a tube from each voxel, return its volume and end surface."""

    domain = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="mask of the calculation domain: one inside, zero outside")
    grad_field = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="scalar field whose gradient is to be advected along")
    divergence = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="divergence of the normalized vector field")
    step_size = Float(
        0.03, output=False, optional=True,
        desc="size of the advection step (millimetres)")
    upfield = Bool(
        False, optional=True,
        desc="Direction of advection (upfield if True, downfield if False)")
    max_dist = Float(
        6, output=False, optional=True,
        desc="maximum advection distance (millimetres)")
    domain_type = Enum(
        "interpolated", values=("boolean", "interpolated"), output=False,
        desc="interpolation type for the domain")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    output_volumes = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output volume containing the tubes' volume")
    output_surfaces = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output volume containing the tubes' end surface")

    def get_commandline(self):
        command_step_size = ((-self.step_size) if self.upfield
                             else self.step_size)
        args = [
            "ylAdvectTubes",
            "--domain", self.domain,
            "--grad-field", self.grad_field,
            "--divergence", self.divergence,
            "--step-size", repr(command_step_size),
            "--max-dist", repr(self.max_dist),
            "--domain-type", self.domain_type,
            "--verbose", str(self.verbosity),
            "--output-volumes", self.output_volumes,
            "--output-surfaces", self.output_surfaces]
        return args


class EuclideanAdvectionAlongGradient(capsul.api.Process):
    """Measure the Euclidean length of an advection path."""

    domain = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        doc="mask of the calculation domain: one inside, zero outside")
    grad_field = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        doc="scalar field whose gradient is to be advected along")
    step_size = Float(
        0.03, output=False, optional=True,
        doc="size of the advection step (millimetres)")
    upfield = Bool(
        False, optional=True,
        desc="Direction of advection (upfield if True, downfield if False)")
    max_dist = Float(
        6, output=False, optional=True,
        doc="maximum advection distance (millimetres)")
    domain_type = Enum(
        "interpolated", values=("boolean", "interpolated"), output=False,
        doc="interpolation type for the domain")
    verbosity = Int(1, output=False, optional=True, doc="Verbosity level")

    output_length = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        doc="output volume containing the length of the advection path")

    def get_commandline(self):
        command_step_size = ((-self.step_size) if self.upfield
                             else self.step_size)
        args = [
            "ylAdvectEuclidean",
            "--domain", self.domain,
            "--grad-field", self.grad_field,
            "--step-size", repr(command_step_size),
            "--max-dist", repr(self.max_dist),
            "--domain-type", self.domain_type,
            "--verbose", str(self.verbosity),
            "--output-length", self.output_length]
        return args


class ImageArithmetic2Inputs(capsul.api.Process):
    """Compute arithmetic from 2 input images"""

    input_image_1 = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image I1")
    input_image_2 = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image I2")
    formula = Str(
        Undefined, output=False,
        desc="arithmetic formula referring to I1 and I2")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="result of the arithmetic"
    )

    def get_commandline(self):
        return [
            "cartoLinearComb.py",
            "-f", self.formula,
            "-i", self.input_image_1,
            "-i", self.input_image_2,
            "-o", self.output_image]


class MergeImagesOneToOne(capsul.api.Process):
    """Merge values into an image using a mask image."""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    mask_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="mask image (must have an integer voxel type)")
    label = Int(
        Undefined, output=False,
        desc="only label of the mask image to take into account")
    value = Float(
        Undefined, output=False,
        desc="replacement value")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output image"
    )

    def get_commandline(self):
        return [
            "AimsMerge",
            "--verbose", "0",
            "-m", "oo",
            "-l", str(self.label),
            "-v", repr(self.value),
            "-i", self.input_image,
            "-M", self.mask_image,
            "-o", self.output_image]

class VolumeSink(capsul.api.Process):
    """Use this process to ignore a mandatory output."""

    file = File(Undefined, allowed_extensions=VOLUME_EXTENSIONS,
                desc="Volume file to be ignored")

    def _run_process(self):
        pass
