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

import enum
import math
import os
import subprocess

import capsul.api
from soma.controller import field, File, Literal, undefined


VOLUME_EXTENSIONS = ['.nii.gz', '.hdr', '.nii', '.ima']


class Laplacian(capsul.api.Process):
    """Solve the Laplacian model in the cortex"""

    classif: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")
    precision: float = field(
        default=0.001,
        doc="target maximum relative error in first-order finite differences")
    typical_cortical_thickness: float = field(
        default=3.0,
        doc="typical thickness of the cortex (mm), used for accelerating "
        "convergence")
    verbosity: int = field(default=1, desc="Verbosity level")

    laplace_field: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output pseudo-temperature field (from 0 in CSF to 1 in the "
        "white matter)"
    )

    def execute(self, context):
        cmd = [
            "ylLaplacian",
            "--classif", self.classif,
            "--output", self.laplace_field,
            "--precision", repr(self.precision),
            "--typical-cortical-thickness",
            repr(self.typical_cortical_thickness),
            "--verbose", str(self.verbosity),
        ]
        subprocess.check_call(cmd)


class IsoCurvature(capsul.api.Process):
    """Compute the curvature of isosurfaces"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image volume (scalar field)")
    # modes mean, geom, pri1, pri2 are currently unimplemented
    mode: Literal["sum"] = field(
        default="sum",
        doc="type of curvature to compute")
    verbosity: int = field(default=1, doc="Verbosity level")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output image volume containing the curvature of the isosurfaces "
        "of the input field"
    )

    def execute(self, context):
        cmd = [
            "ylIsoCurvature",
            "--input", self.input_image,
            "--mode", self.mode,
            "--output", self.output_image,
            "--verbose", str(self.verbosity),
        ]
        subprocess.check_call(cmd)


class RemoveNaN(capsul.api.Process):
    """Remove NaN values from an image"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    value: float = field(
        default=0.0,
        doc="replacement value")
    percentage: bool = field(
        default=True,
        doc="interpret value as a percentage of the image intensity range")

    output_image: File = field(
        write=True,
        extensions=VOLUME_EXTENSIONS,
        doc="output image"
    )

    def execute(self, context):
        cmd = [
            "AimsRemoveNaN",
            "--verbose", "0",
            "-i", self.input_image,
            "-np", str(self.percentage),
            "--value", repr(self.value),
            "-o", self.output_image
        ]
        subprocess.check_call(cmd)


class MedianFilter(capsul.api.Process):
    """Median filter smoothing"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    x_size: int = field(default=3,
                        doc="X size of the filter mask")
    y_size: int = field(default=3,
                        doc="Y size of the filter mask")
    z_size: int = field(default=3,
                        doc="Z size of the filter mask")

    output_image: File = field(
        write=True,
        extensions=VOLUME_EXTENSIONS,
        doc="median-filtered image"
    )

    def execute(self, context):
        cmd = [
            "AimsMedianSmoothing",
            "--verbose", "0",
            "--input", self.input_image,
            "--dx", str(self.x_size),
            "--dy", str(self.y_size),
            "--dz", str(self.z_size),
            "--output", self.output_image
        ]
        subprocess.check_call(cmd)


class GaussianSmoothing(capsul.api.Process):
    """3D Gaussian smoothing filter using the recursive Deriche method"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    xsigma: float = field(
        default=None,
        doc="X standard deviation of the gaussian filter "
        "[default=largest voxel size]")
    ysigma: float = field(
        default=None,
        doc="Y standard deviation of the gaussian filter "
        "[default=largest voxel size]")
    zsigma: float = field(
        default=None,
        doc="Z standard deviation of the gaussian filter "
        "[default=largest voxel size]")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="Gaussian-filtered image"
    )

    def execute(self, context):
        sigma_args = []
        if self.xsigma is not None:
            sigma_args += ["--xsigma", str(self.xsigma)]
        if self.ysigma is not None:
            sigma_args += ["--ysigma", str(self.ysigma)]
        if self.zsigma is not None:
            sigma_args += ["--zsigma", str(self.zsigma)]
        cmd = [
            "AimsGaussianSmoothing",
            "--input", self.input_image
        ] + sigma_args + [
            "--output", self.output_image
        ]
        subprocess.check_call(cmd)


class BinarizeCortex(capsul.api.Process):
    """Extract a binary image (0/1) of the cortex"""

    classif: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="binary image of the cortex (1 in the cortex, 0 elsewhere)")

    def execute(self, context):
        cmd = [
            "AimsThreshold",
            "--verbose", "0",
            "-b",
            "--fg", "1",
            "-m", "eq",
            "-t", "100",
            "--input", self.classif,
            "--output", self.output_image,
        ]
        subprocess.check_call(cmd)


class DomainTypeEnum(str, enum.Enum):
    interpolated = "interpolated"
    boolean = "boolean"


class AdvectTubesAlongGradient(capsul.api.Process):
    """Advect a tube from each voxel, return its volume and end surface."""

    domain: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="mask of the calculation domain: one inside, zero outside")
    grad_field: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="scalar field whose gradient is to be advected along")
    divergence: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="divergence of the normalized vector field")
    step_size: float = field(
        default=0.03,
        doc="size of the advection step (millimetres)")
    upfield: bool = field(
        doc="Direction of advection (upfield if True, downfield if False)")
    max_dist: float = field(
        default=6,
        doc="maximum advection distance (millimetres)")
    domain_type: DomainTypeEnum = field(
        default="interpolated",
        doc="interpolation type for the domain")
    verbosity: int = field(default=1, doc="Verbosity level")

    output_volumes: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output volume containing the tubes' volume")
    output_surfaces: File = field(
        write=True, extensions=VOLUME_EXTENSIONS, optional=True,
        doc="output volume containing the tubes' end surface")

    def execute(self, context):
        command_step_size = ((-self.step_size) if self.upfield
                             else self.step_size)
        cmd = [
            "ylAdvectTubes",
            "--domain", self.domain,
            "--grad-field", self.grad_field,
            "--divergence", self.divergence,
            "--step-size", repr(command_step_size),
            "--max-dist", repr(self.max_dist),
            "--domain-type", self.domain_type.value,
            "--verbose", str(self.verbosity),
            "--output-volumes", self.output_volumes,
        ]
        if self.output_surfaces is not undefined:
            cmd += ["--output-surfaces", self.output_surfaces]

        subprocess.check_call(cmd)


class EuclideanAdvectionAlongGradient(capsul.api.Process):
    """Measure the Euclidean length of an advection path."""

    domain: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="mask of the calculation domain: one inside, zero outside")
    grad_field: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="scalar field whose gradient is to be advected along")
    step_size: float = field(
        default=0.03,
        doc="size of the advection step (millimetres)")
    upfield: bool = field(
        doc="Direction of advection (upfield if True, downfield if False)")
    max_dist: float = field(
        default=6,
        doc="maximum advection distance (millimetres)")
    domain_type: DomainTypeEnum = field(
        default="interpolated",
        doc="interpolation type for the domain")
    verbosity: int = field(default=1, doc="Verbosity level")

    output_length: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output volume containing the length of the advection path")

    def execute(self, context):
        command_step_size = ((-self.step_size) if self.upfield
                             else self.step_size)
        cmd = [
            "ylAdvectEuclidean",
            "--domain", self.domain,
            "--grad-field", self.grad_field,
            "--step-size", repr(command_step_size),
            "--max-dist", repr(self.max_dist),
            "--domain-type", self.domain_type.value,
            "--verbose", str(self.verbosity),
            "--output-length", self.output_length
        ]
        subprocess.check_call(cmd)


class PostProcessEquivolumetricDepth(capsul.api.Process):
    """Post-process an equivolumetric depth image.

    - Set the outside of the brain (CSF) to 0.0
    - Set the white matter to 1.0
    - Set various Nifti header fields
    """

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image of equivolumetric depth")
    classif: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output image"
    )

    def execute(self, context):
        cmd = [
            "ylPostProcessEquivolumetricDepth",
            self.input_image,
            self.classif,
            self.output_image
        ]
        subprocess.check_call(cmd)


class ImageArithmetic2Inputs(capsul.api.Process):
    """Compute arithmetic from 2 input images"""

    input_image_1: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image I1")
    input_image_2: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image I2")
    formula: str = field(
        doc="arithmetic formula referring to I1 and I2")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="result of the arithmetic"
    )

    def execute(self, context):
        cmd = [
            "cartoLinearComb.py",
            "-f", self.formula,
            "-i", self.input_image_1,
            "-i", self.input_image_2,
            "-o", self.output_image
        ]
        subprocess.check_call(cmd)


class MergeImagesOneToOne(capsul.api.Process):
    """Merge values into an image using a mask image."""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    mask_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="mask image (must have an integer voxel type)")
    label_to_replace: int = field(
        doc="only label of the mask image to take into account")
    value: float = field(
        doc="replacement value")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output image"
    )

    def execute(self, context):
        cmd = [
            "AimsMerge",
            "--verbose", "0",
            "-m", "oo",
            "-l", str(self.label_to_replace),
            "-v", repr(self.value),
            "-i", self.input_image,
            "-M", self.mask_image,
            "-o", self.output_image
        ]
        subprocess.check_call(cmd)


class EuclideanUpwindingAlongGradient(capsul.api.Process):
    """Compute distance to a boundary along the gradient of a scalar field."""

    domain: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="label image defining the computation domain")
    scalar_field: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="scalar field whose gradient is used as the integration "
        "direction")
    downfield: bool = field(
        doc="work on inverted field (downfield instead of upfield)")
    domain_label: int = field(
        default=100,
        doc="label of the propagation domain")
    origin_label: int = field(
        default=0,
        doc="label of the origin object")
    verbosity: int = field(default=1, doc="Verbosity level")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output volume containing the distance")

    def execute(self, context):
        cmd = [
            "ylUpwindDistance",
            "--domain", self.domain,
            "--field", self.scalar_field,
            "--invert", str(self.downfield),
            "--domain-label", str(self.domain_label),
            "--origin-label", str(self.origin_label),
            "--verbose", str(self.verbosity),
            "--output", self.output
        ]
        subprocess.check_call(cmd)


class Distmaps(capsul.api.Process):
    """Compute distance maps to the boundaries of the cortex"""

    classif: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")

    distwhite: File = field(
        default=os.devnull,
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="signed Euclidean distance to the white matter interface"
    )
    distCSF: File = field(
        default=os.devnull,
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="signed Euclidean distance to the CSF interface"
    )
    classif_with_outer_boundaries: File = field(
        default=os.devnull,
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex with labelled boundaries "
        "(50 on the CSF, 150 on the white matter)")

    def execute(self, context):
        cmd = [
            "ylDistmaps",
            self.classif,
            self.distwhite,
            self.distCSF,
            self.classif_with_outer_boundaries
        ]
        subprocess.check_call(cmd)


class ThresholdModeEnum(str, enum.Enum):
    eq = "eq"
    di = "di"
    lt = "lt"
    le = "le"
    gt = "gt"
    ge = "ge"


class ImageSingleThreshold(capsul.api.Process):
    """Threshold an image"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    binary: bool = field(
        default=False,
        doc="return a binary result as int16")
    fg: int = field(
        default=32767,
        doc="foreground value set on thresholded in voxels in binary mode")
    threshold: float = field(
        doc="value of the threshold")
    mode: ThresholdModeEnum = field(
        doc="""\
thresholding type
    lt   --> lower than
    le   --> lower or equal to
    gt   --> greater than
    ge   --> greater or equal to
    eq   --> equal to
    di   --> differ
""")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="thresholded image")

    def execute(self, context):
        cmd = [
            "AimsThreshold",
            "--verbose", "0",
            "-b", str(self.binary),
            "-m", self.mode.value,
            "-t", repr(self.threshold),
            "--input", self.input_image,
            "--output", self.output_image
        ]
        if self.binary:
            cmd += ["--fg", str(self.fg)]
        subprocess.check_call(cmd)


class LabelEachVoxel(capsul.api.Process):
    """Assign a unique label to each voxel of a mask"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input mask")
    first_label: int = field(
        default=1,
        doc="assign labels starting with this value")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output label volume with S32 datatype")

    def execute(self, context):
        cmd = [
            "ylLabelEachVoxel",
            "--first-label", str(self.first_label),
            "--input", self.input_image,
            "--output", self.output_image
        ]
        subprocess.check_call(cmd)


class VolumeDataTypeEnum(str, enum.Enum):
    CDOUBLE = "CDOUBLE"
    CFLOAT = "CFLOAT"
    DOUBLE = "DOUBLE"
    FLOAT = "FLOAT"
    HSV = "HSV"
    POINT3DF = "POINT3DF"
    RGB = "RGB"
    RGBA = "RGBA"
    S16 = "S16"
    S32 = "S32"
    S8 = "S8"
    U16 = "U16"
    U32 = "U32"
    U8 = "U8"
    VECTOR_OF_3_SHORT = "VECTOR_OF_3_SHORT"
    VECTOR_OF_6_FLOAT = "VECTOR_OF_6_FLOAT"


class ConvertDataType(capsul.api.Process):
    """Convert the data type of an image"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    data_type: VolumeDataTypeEnum = field(
        doc="output data type")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output label volume with S32 datatype")

    def execute(self, context):
        cmd = [
            "AimsFileConvert",
            "--type", self.data_type.value,
            "--input", self.input_image,
            "--output", self.output_image
        ]
        subprocess.check_call(cmd)


class MergeImagesAllToOne(capsul.api.Process):
    """Merge values into an image using a mask image."""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    mask_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="mask image (must have an integer voxel type)")
    value: float = field(
        doc="replacement value")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output image"
    )

    def execute(self, context):
        cmd = [
            "AimsMerge",
            "--mode", "ao",
            "--value", repr(self.value),
            "--input", self.input_image,
            "--Mask", self.mask_image,
            "--output", self.output_image
        ]
        subprocess.check_call(cmd)


class MergeImagesSameValues(capsul.api.Process):
    """Merge values into an image using a mask image."""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input image")
    mask_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="mask image (must have an integer voxel type)")

    output_image: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output image"
    )

    def execute(self, context):
        cmd = [
            "AimsMerge",
            "--mode", "sv",
            "--input", self.input_image,
            "--Mask", self.mask_image,
            "--output", self.output_image
        ]
        subprocess.check_call(cmd)


class PropagateAlongFieldGradient(capsul.api.Process):
    """Propagate labels along the gradient of a scalar field."""

    seeds: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="""\
volume of labels (either S16 or S32):
    - positive labels are seeds,
    - zero is the region of propagation,
    - negative labels are forbidden regions.
""")
    target_label: int = field(
        default=0,
        doc="voxels having this label are used as advection starting points")
    grad_field: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="scalar field whose gradient is to be advected along")
    step_size: float = field(
        default=0.03,
        doc="size of the advection step (millimetres)")
    upfield: bool = field(
        doc="Direction of advection (upfield if True, downfield if False)")
    max_dist: float = field(
        default=6,
        doc="maximum advection distance (millimetres)")
    verbosity: int = field(default=1, doc="Verbosity level")

    output_labels: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output the propagated labels")
    dest_points: File = field(
        write=True, optional=True,
        extensions=VOLUME_EXTENSIONS,
        doc="output the destination points for each propagated voxel")

    def execute(self, context):
        command_step_size = ((-self.step_size) if self.upfield
                             else self.step_size)
        cmd = [
            "ylPropagateAlongField",
            "--seeds", self.seeds,
            "--grad-field", self.grad_field,
            "--target-label", str(self.target_label),
            "--step", repr(command_step_size),
            "--max-iter", str(int(math.ceil(self.max_dist / self.step_size))),
            "--verbose", str(self.verbosity),
            "--output", self.output_labels,
        ]
        if self.dest_points is not undefined:
            cmd += ["--dest-points", self.dest_points]
        subprocess.check_call(cmd)


class GetExchangedPropagationVolume(capsul.api.Process):
    """Get a volume of exchanged propagation labels"""

    classif_with_outer_boundaries: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter, 50 on the CSF border, 150 on the white matter "
        "border)")
    CSF_labels_on_white: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="labels of the CSF projected onto the white matter boundary")
    white_labels_on_CSF: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="labels of the white matter projected onto the CSF boundary")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="volume where each interface is labelled with connected "
        "components facing the same voxels of the other interface")

    def execute(self, context):
        cmd = [
            "ylGetExchangedPropvol",
            self.classif_with_outer_boundaries,
            self.CSF_labels_on_white,
            self.white_labels_on_CSF,
            self.output
        ]
        subprocess.check_call(cmd)


class RelabelConjunction(capsul.api.Process):
    """Assign new labels to voxels that have the same pair of labels"""

    labels1: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input label image")
    labels2: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input label image")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output label image")

    def execute(self, context):
        cmd = [
            "ylRelabelConjunction",
            self.labels1,
            self.labels2,
            self.output
        ]
        subprocess.check_call(cmd)


# TODO restore enum
# class ConnectivityTypeEnum(str, enum.Enum):
#     c26: "26"
#     c4xy: "4xy"  # noqa: F722
#     c4xz: "4xz"  # noqa: F722
#     c4yz: "4yz"  # noqa: F722
#     c6: "6"
#     c8xy: "8xy"  # noqa: F722
#     c8xz: "8xz"  # noqa: F722
#     c8yz: "8yz"  # noqa: F722
#     c18: "18"
ConnectivityTypeEnum = str


class ConnectedComponents(capsul.api.Process):
    """Extract connected components of a labelled volume"""

    input_image: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input label image")
    connectivity: ConnectivityTypeEnum = field(
        default="26",
        doc="connectivity")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output labelled connected components volume")

    def execute(self, context):
        cmd = [
            "AimsConnectComp",
            "--input", self.input_image,
            "--output", self.output,
            "--connectivity", self.connectivity
        ]
        subprocess.check_call(cmd)


class MergeCortexColumnRegions(capsul.api.Process):
    """Aggregate over-segmented cortical traverses."""

    input_traverses: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input label volume")
    proj_csf: File = field(
        optional=True,
        extensions=VOLUME_EXTENSIONS,
        doc="projected coordinates of the CSF surface")
    proj_white: File = field(
        optional=True,
        extensions=VOLUME_EXTENSIONS,
        doc="projected coordinates of the white surface")
    classif: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="classification image of the cortex (100 inside, 0 in CSF, "
        "200 in white matter)")
    goal_diameter: float = field(
        default=0.5,
        doc="goal region diameter (millimetres)")
    verbosity: int = field(default=2, doc="Verbosity level")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output label volume")

    def execute(self, context):
        cmd = [
            "ylMergeCortexColumnRegions",
            "--input", self.input_traverses,
            "--proj-csf", self.proj_csf,
            "--proj-white", self.proj_white,
            "--classif", self.classif,
            "--goal-diameter", repr(self.goal_diameter),
            "--verbose", str(self.verbosity),
            "--output", self.output,
        ]
        subprocess.check_call(cmd)


class Relabel(capsul.api.Process):
    """Assign new consecutive labels to an existing label image"""

    input: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input label image")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output label image")

    def execute(self, context):
        cmd = [
            "ylRelabel",
            self.input,
            self.output
        ]
        subprocess.check_call(cmd)


class RandomizeLabels(capsul.api.Process):
    """Randomize the labels of an image with consecutive labels"""

    input: File = field(
        extensions=VOLUME_EXTENSIONS,
        doc="input label image")

    output: File = field(
        write=True, extensions=VOLUME_EXTENSIONS,
        doc="output label image")

    def execute(self, context):
        cmd = [
            "ylRandomizeLabels",
            self.input,
            self.output
        ]
        subprocess.check_call(cmd)


# TODO delete
if __name__ == '__main__':
    import sys
    sys.stdout.flush()
    from soma.qt_gui.qt_backend import QtGui
    from capsul.qt_gui.widgets import PipelineDeveloperView
    pipeline = capsul.api.Capsul.executable(
        'highres_cortex.capsul.isovolume')
    app = QtGui.QApplication.instance()
    if not app:
        app = QtGui.QApplication(sys.argv)
    view1 = PipelineDeveloperView(pipeline, show_sub_pipelines=True)
    view1.show()
    app.exec_()
