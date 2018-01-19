# -*- coding: utf-8 -*-
# Copyright Forschungszentrum Jülich GmbH (2017, 2018).
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

import math

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
        "sum", output=False, optional=True,
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

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output image"
    )

    def get_commandline(self):
        return [
            "AimsRemoveNaN",
            "--verbose", "0",
            "-i", self.input_image,
            "-np", str(self.percentage),
            "--value", repr(self.value),
            "-o", self.output_image]


class MedianFilter(capsul.api.Process):
    """Median filter smoothing"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    x_size = Int(3, output=False, optional=True,
                 desc="X size of the filter mask")
    y_size = Int(3, output=False, optional=True,
                 desc="Y size of the filter mask")
    z_size = Int(3, output=False, optional=True,
                 desc="Z size of the filter mask")

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
        False, optional=False,
        desc="Direction of advection (upfield if True, downfield if False)")
    max_dist = Float(
        6, output=False, optional=True,
        desc="maximum advection distance (millimetres)")
    domain_type = Enum(
        "interpolated", "boolean", output=False, optional=True,
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
        desc="mask of the calculation domain: one inside, zero outside")
    grad_field = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="scalar field whose gradient is to be advected along")
    step_size = Float(
        0.03, output=False, optional=True,
        desc="size of the advection step (millimetres)")
    upfield = Bool(
        False, optional=False,
        desc="Direction of advection (upfield if True, downfield if False)")
    max_dist = Float(
        6, output=False, optional=True,
        desc="maximum advection distance (millimetres)")
    domain_type = Enum(
        "interpolated", "boolean", output=False, optional=True,
        desc="interpolation type for the domain")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    output_length = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output volume containing the length of the advection path")

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
        # bv_env automatically launches the command through Python on Windows
        return [
            "bv_env",
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


class EuclideanUpwindingAlongGradient(capsul.api.Process):
    """Compute distance to a boundary along the gradient of a scalar field."""

    domain = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="label image defining the computation domain")
    field = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="scalar field whose gradient is used as the integration "
        "direction")
    downfield = Bool(
        False, optional=False,
        desc="work on inverted field (downfield instead of upfield)")
    domain_label = Int(
        100, optional=True,
        desc="label of the propagation domain")
    origin_label = Int(
        0, optional=True,
        desc="label of the origin object")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output volume containing the distance")

    def get_commandline(self):
        return [
            "ylUpwindDistance",
            "--domain", self.domain,
            "--field", self.field,
            "--invert", str(self.downfield),
            "--domain-label", str(self.domain_label),
            "--origin-label", str(self.origin_label),
            "--verbose", str(self.verbosity),
            "--output", self.output]


class Distmaps(capsul.api.Process):
    """Compute distance maps to the boundaries of the cortex"""

    classif = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")

    distwhite = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="signed Euclidean distance to the white matter interface"
    )
    distCSF = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="signed Euclidean distance to the CSF interface"
    )
    classif_with_outer_boundaries = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="classification image of the cortex with labelled boundaries "
        "(50 on the CSF, 150 on the white matter)")

    def get_commandline(self):
        # bv_env automatically launches the command through Python on Windows
        return [
            "bv_env",
            "ylDistmaps",
            self.classif,
            self.distwhite,
            self.distCSF,
            self.classif_with_outer_boundaries
        ]


class ImageSingleThreshold(capsul.api.Process):
    """Threshold an image"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    binary = Bool(
        False, output=False, optional=True,
        desc="return a binary result as int16")
    fg = Int(
        32767, output=False, optional=True,
        desc="foreground value set on thresholded in voxels in binary mode")
    threshold = Float(
        Undefined, output=False,
        desc="value of the threshold")
    mode = Enum(
        "eq", "di", "lt", "le", "gt", "ge", output=False,
        desc="""\
thresholding type
    lt   --> lower than
    le   --> lower or equal to
    gt   --> greater than
    ge   --> greater or equal to
    eq   --> equal to
    di   --> differ
""")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="thresholded image")

    def get_commandline(self):
        cmd = [
            "AimsThreshold",
            "--verbose", "0",
            "-b", str(self.binary),
            "-m", self.mode,
            "-t", repr(self.threshold),
            "--input", self.input_image,
            "--output", self.output_image
        ]
        if self.binary:
            cmd += ["--fg", str(self.fg)]
        return cmd


class LabelEachVoxel(capsul.api.Process):
    """Assign a unique label to each voxel of a mask"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input mask")
    first_label = Int(
        1, output=False, optional=True,
        desc="assign labels starting with this value")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output label volume with S32 datatype")

    def get_commandline(self):
        return [
            "ylLabelEachVoxel",
            "--first-label", str(self.first_label),
            "--input", self.input_image,
            "--output", self.output_image
        ]


class ConvertDataType(capsul.api.Process):
    """Convert the data type of an image"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    data_type = Enum(
        "CDOUBLE", "CFLOAT", "DOUBLE", "FLOAT", "HSV", "POINT3DF", "RGB",
        "RGBA", "S16", "S32", "S8", "U16", "U32", "U8", "VECTOR_OF_3_SHORT",
        "VECTOR_OF_6_FLOAT", output=False,
        desc="output data type")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output label volume with S32 datatype")

    def get_commandline(self):
        return [
            "AimsFileConvert",
            "--type", self.data_type,
            "--input", self.input_image,
            "--output", self.output_image
        ]


class MergeImagesAllToOne(capsul.api.Process):
    """Merge values into an image using a mask image."""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    mask_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="mask image (must have an integer voxel type)")
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
            "--mode", "ao",
            "--value", repr(self.value),
            "--input", self.input_image,
            "--Mask", self.mask_image,
            "--output", self.output_image
        ]


class MergeImagesSameValues(capsul.api.Process):
    """Merge values into an image using a mask image."""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input image")
    mask_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="mask image (must have an integer voxel type)")

    output_image = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output image"
    )

    def get_commandline(self):
        return [
            "AimsMerge",
            "--mode", "sv",
            "--input", self.input_image,
            "--Mask", self.mask_image,
            "--output", self.output_image
        ]


class PropagateAlongFieldGradient(capsul.api.Process):
    """Propagate labels along the gradient of a scalar field."""

    seeds = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="""\
volume of labels (either S16 or S32):
    - positive labels are seeds,
    - zero is the region of propagation,
    - negative labels are forbidden regions.
""")
    target_label = Int(
        0, output=False, optional=True,
        desc="voxels having this label are used as advection starting points")
    grad_field = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="scalar field whose gradient is to be advected along")
    step_size = Float(
        0.03, output=False, optional=True,
        desc="size of the advection step (millimetres)")
    upfield = Bool(
        False, optional=False,
        desc="Direction of advection (upfield if True, downfield if False)")
    max_dist = Float(
        6, output=False, optional=True,
        desc="maximum advection distance (millimetres)")
    verbosity = Int(1, output=False, optional=True, desc="Verbosity level")

    output_labels = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output the propagated labels")
    dest_points = File(
        Undefined, output=True, optional=True,
        allowed_extensions=VOLUME_EXTENSIONS,
        desc="output the destination points for each propagated voxel")

    def get_commandline(self):
        command_step_size = ((-self.step_size) if self.upfield
                             else self.step_size)
        args = [
            "ylPropagateAlongField",
            "--seeds", self.seeds,
            "--grad-field", self.grad_field,
            "--target-label", str(self.target_label),
            "--step", repr(command_step_size),
            "--max-iter", str(int(math.ceil(self.max_dist / self.step_size))),
            "--verbose", str(self.verbosity),
            "--output", self.output_labels,
        ]
        if self.dest_points is not Undefined:
            args += ["--dest-points", self.dest_points]
        return args


class GetExchangedPropagationVolume(capsul.api.Process):
    """Get a volume of exchanged propagation labels"""

    classif_with_outer_boundaries = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter, 50 on the CSF border, 150 on the white matter "
        "border)")
    CSF_labels_on_white = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="labels of the CSF projected onto the white matter boundary")
    white_labels_on_CSF = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="labels of the white matter projected onto the CSF boundary")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="volume where each interface is labelled with connected "
        "components facing the same voxels of the other interface")

    def get_commandline(self):
        # bv_env automatically launches the command through Python on Windows
        return [
            "bv_env",
            "ylGetExchangedPropvol",
            self.classif_with_outer_boundaries,
            self.CSF_labels_on_white,
            self.white_labels_on_CSF,
            self.output
        ]


class RelabelConjunction(capsul.api.Process):
    """Assign new labels to voxels that have the same pair of labels"""

    labels1 = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input label image")
    labels2 = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input label image")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output label image")

    def get_commandline(self):
        # bv_env automatically launches the command through Python on Windows
        return [
            "bv_env",
            "ylRelabelConjunction",
            self.labels1,
            self.labels2,
            self.output
        ]


class ConnectedComponents(capsul.api.Process):
    """Extract connected components of a labelled volume"""

    input_image = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input label image")
    connectivity = Enum(
        "26", "4xy", "4xz", "4yz", "6", "8xy", "8xz", "8yz", "18",
        output=False, optional=True,
        desc="connectivity")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output labelled connected components volume")

    def get_commandline(self):
        return [
            "AimsConnectComp",
            "--input", self.input_image,
            "--output", self.output,
            "--connectivity", self.connectivity,
        ]


class MergeCortexColumnRegions(capsul.api.Process):
    """Aggregate over-segmented cortical traverses."""

    input_traverses = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input label volume")
    proj_csf = File(
        Undefined, output=False, optional=True,
        allowed_extensions=VOLUME_EXTENSIONS,
        desc="projected coordinates of the CSF surface")
    proj_white = File(
        Undefined, output=False, optional=True,
        allowed_extensions=VOLUME_EXTENSIONS,
        desc="projected coordinates of the white surface")
    classif = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")
    goal_diameter = Float(
        0.5, output=False, optional=True,
        desc="goal region diameter (millimetres)")
    verbosity = Int(2, output=False, optional=True, desc="Verbosity level")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output label volume")

    def get_commandline(self):
        args = [
            "ylMergeCortexColumnRegions",
            "--input", self.input_traverses,
            "--proj-csf", self.proj_csf,
            "--proj-white", self.proj_white,
            "--classif", self.classif,
            "--goal-diameter", repr(self.goal_diameter),
            "--verbose", str(self.verbosity),
            "--output", self.output,
        ]
        return args


class Relabel(capsul.api.Process):
    """Assign new consecutive labels to an existing label image"""

    input = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input label image")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output label image")

    def get_commandline(self):
        # bv_env automatically launches the command through Python on Windows
        return [
            "bv_env",
            "ylRelabel",
            self.input,
            self.output
        ]


class RandomizeLabels(capsul.api.Process):
    """Randomize the labels of an image with consecutive labels"""

    input = File(
        Undefined, output=False, allowed_extensions=VOLUME_EXTENSIONS,
        desc="input label image")

    output = File(
        Undefined, output=True, allowed_extensions=VOLUME_EXTENSIONS,
        desc="output label image")

    def get_commandline(self):
        # bv_env automatically launches the command through Python on Windows
        return [
            "bv_env",
            "ylRandomizeLabels",
            self.input,
            self.output
        ]
