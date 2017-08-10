#! /bin/sh -e
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

AimsThreshold -b -m eq -t 50 \
    -i ../classif_with_outer_boundaries.nii.gz \
    -o CSF_interface.nii
AimsThreshold -b -m eq -t 150 \
    -i ../classif_with_outer_boundaries.nii.gz \
    -o white_interface.nii
ylLabelEachVoxel --verbose \
    -i CSF_interface.nii.gz \
    -o CSF_labelled_interface.nii \
    --first-label 100000001
ylLabelEachVoxel --verbose \
    -i white_interface.nii.gz \
    -o white_labelled_interface.nii \
    --first-label 200000001

AimsThreshold -b --fg -1 -m di -t 100 \
    -i ../classif.nii.gz \
    -o negative_outside_cortex.nii
AimsFileConvert -t S32 \
    -i negative_outside_cortex.nii \
    -o negative_outside_cortex_S32.nii

AimsMerge -m sv \
    -i negative_outside_cortex_S32.nii \
    -M CSF_labelled_interface.nii \
    -o CSF_labelled_interface_negative_outside.nii
AimsMerge -m ao -v 200000000 \
    -i CSF_labelled_interface_negative_outside.nii \
    -M white_labelled_interface.nii \
    -o propvol_CSF_labels.nii.gz

AimsMerge -m sv \
    -i negative_outside_cortex_S32.nii \
    -M white_labelled_interface.nii \
    -o white_labelled_interface_negative_outside.nii
AimsMerge -m ao -v 100000000 \
    -i white_labelled_interface_negative_outside.nii \
    -M CSF_labelled_interface.nii \
    -o propvol_white_labels.nii.gz


ylPropagateAlongField --verbose \
    --grad-field ../heat/heat.nii.gz \
    --seeds propvol_CSF_labels.nii.gz \
    --step -0.05 \
    --target-label 200000000 \
    --output heat_CSF_labels_on_white.nii.gz
ylPropagateAlongField --verbose \
    --grad-field ../heat/heat.nii.gz \
    --seeds propvol_white_labels.nii.gz \
    --step 0.05 \
    --target-label 100000000 \
    --output heat_white_labels_on_CSF.nii.gz

python "$(dirname "$0")"/get_exchanged_propvol.py  # -> exchanged_propvol.nii.gz

# Why is the previous step necessary?
#
# The obvious alternative is to do exactly as described in the OHBM paper: do
# the projections on the original labels of each voxel.
#
# The previous case aggregates the adjacent voxels of one interface that point
# towards the same voxel on the other interface. This reduces
# over-segmentation.
#
# Another way of reducing over-segmentation would be to aggregate together
# voxels that have one projection in common, instead of both (see conjunction
# step later on). But this introduces the problem of transitivity. This was
# investigated previously on the ferret data (under the name Billiard), but was
# considered a dead-end and the above solution seems to solve this problem most
# efficiently.


# There is a problem with the propagation of labels: the step size is fixed,
# which means that sometimes the point can skip the corner of a voxel, and thus
# go directly from a bulk voxel to an outside voxel. In this case it is
# recorded as a "dead-end" advection path, no resulting label is recorded and
# it appears as zero in the result.
#
# This problem also appears in the previous "exchange" step, but is mitigated
# by the subsequent connex component detection (each failed propagation is
# assigned a different label).
#
# Quick fix: fix the conjunction step to not aggregate zeros.
#
# TODO: the proper way to fix this would be to force the advection path to
# respect the boundaries of voxels, so that the corner of voxels cannot be
# skipped over. This would also prevent the advection path from crossing the
# thin CSF surface within the sulcus (comes from skeleton).

# I could take into account the fake cortex–CSF interface that exists at the
# cut plane, by assigning it a special label (e.g. 500000000) in the
# exchanged_propvol label. It would then need to be treated specially: any
# voxel that projects onto this label would be excluded from the region list,
# and thus would not take part in the merging step. This would prevent the
# creation of regions that connect to this spurious surface, but this would not
# prevent the nearby regions from being deformed by the perturbation of the
# field. It would thus probably be overkill to implement this special case.
# Care is needed when dealing with regions close to the cut plane anyway.

AimsMerge -m oo -l 150 -v 0 \
    -i exchanged_propvol.nii.gz \
    -M ../classif_with_outer_boundaries.nii.gz \
    -o ./exchanged_labels_on_CSF.nii
AimsMerge -m oo -l 50 -v 0 \
    -i ./exchanged_propvol.nii.gz \
    -M ../classif_with_outer_boundaries.nii.gz \
    -o ./exchanged_labels_on_white.nii

ylPropagateAlongField --verbose \
    --grad-field ../heat/heat.nii.gz \
    --seeds exchanged_labels_on_CSF.nii \
    --step -0.05 \
    --target-label 0 \
    --output heat_CSF_on_bulk.nii.gz \
    --dest-points heat_CSF_points_on_bulk.nii.gz
ylPropagateAlongField --verbose \
    --grad-field ../heat/heat.nii.gz \
    --seeds exchanged_labels_on_white.nii \
    --step 0.05 \
    --target-label 0 \
    --output heat_white_on_bulk.nii.gz \
    --dest-points heat_white_points_on_bulk.nii.gz

python "$(dirname "$0")"/relabel_conjunction.py  # -> ./conjunction.nii.gz

AimsConnectComp -c 26 \
    -i conjunction.nii.gz \
    -o conjunction_connected.nii.gz

ylMergeCortexColumnRegions --verbose 2 \
    -i conjunction_connected.nii.gz \
    -o merged.nii \
    --proj-csf heat_CSF_points_on_bulk.nii.gz \
    --proj-white heat_white_points_on_bulk.nii.gz \
    --classif ../classif.nii.gz \
    --goal-diameter 1
python "$(dirname "$0")"/relabel.py
python "$(dirname "$0")"/randomize_labels.py
