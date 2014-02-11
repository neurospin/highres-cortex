#! /bin/sh -e

AimsThreshold -b -m eq -t 0 \
    -i ../dist/distCSF.nii.gz \
    -o CSF_interface.nii
AimsThreshold -b -m eq -t 0 \
    -i ../dist/distwhite.nii.gz \
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
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --seeds propvol_CSF_labels.nii.gz \
    --step -0.05 \
    --target-label 200000000 \
    --output heat_CSF_labels_on_white.nii.gz
ylPropagateAlongField --verbose \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --seeds propvol_white_labels.nii.gz \
    --step 0.05 \
    --target-label 100000000 \
    --output heat_white_labels_on_CSF.nii.gz

python get_exchanged_propvol.py  # -> exchanged_propvol.nii.gz

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

ylPropagateAlongField --verbose \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --seeds exchanged_propvol.nii.gz \
    --step -0.05 \
    --target-label 0 \
    --output heat_CSF_on_bulk_raw.nii \
    --dest-points heat_CSF_points_on_bulk.nii.gz
ylPropagateAlongField --verbose \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --seeds exchanged_propvol.nii.gz \
    --step 0.05 \
    --target-label 0 \
    --output heat_white_on_bulk_raw.nii \
    --dest-points heat_white_points_on_bulk.nii.gz

AimsThreshold -b -m eq -t 100 \
    -i ../classif.nii.gz \
    -o cortex_mask.nii

AimsMask \
    -i heat_CSF_on_bulk_raw.nii \
    -m cortex_mask.nii \
    -o heat_CSF_on_bulk.nii.gz
AimsMask \
    -i heat_white_on_bulk_raw.nii \
    -m cortex_mask.nii \
    -o heat_white_on_bulk.nii.gz

python relabel_conjunction.py  # -> ./conjunction.nii.gz

#AimsConnectComp -i conjunction.nii.gz -o conjunction_connected.nii.gz

ylMergeCortexColumnRegions --verbose 2 \
    -i conjunction.nii.gz \
    -o merged.nii.gz \
    --proj-csf heat_CSF_points_on_bulk.nii.gz \
    --proj-white ./heat_white_points_on_bulk.nii.gz \
    --goal-diametre 4 \
    --max-thickness 6
python relabel.py
python randomize_labels.py


N=29719 ; python \
    /volatile/yl232319/bv/src/perso/leprince/tools/proj_point_cloud.py \
    heat_CSF_points_on_bulk.nii.gz \
    merged_randomized.nii.gz \
    csf_projected_points.gii \
    $N \
&& python \
    /volatile/yl232319/bv/src/perso/leprince/tools/proj_point_cloud.py \
    heat_white_points_on_bulk.nii.gz \
    merged_randomized.nii.gz \
    white_projected_points.gii \
    $N \
&& AimsThreshold -b -m eq -t $N \
    -i merged_randomized.nii.gz \
    -o focused_region.nii \
&& AimsGraphConvert --bucket \
    -i focused_region.nii \
    -o focused_region.arg
