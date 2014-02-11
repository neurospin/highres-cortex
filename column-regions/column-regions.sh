AimsThreshold -b -i ../dist/distmap_from_CSF.nii.gz -m eq -t 0 -o "$TMP"/CSF_interface.nii.gz
AimsThreshold -b -i ../dist/distmap_from_white.nii.gz -m eq -t 0 -o "$TMP"/white_interface.nii.gz
ylLabelEachVoxel --verbose -i "$TMP"/CSF_interface.nii.gz -o "$TMP"/CSF_labelled_interface.nii.gz --first-label 100000001
ylLabelEachVoxel --verbose -i "$TMP"/white_interface.nii.gz -o "$TMP"/white_labelled_interface.nii.gz --first-label 200000001

AimsReplaceLevel -i ../classif.nii.gz -o "$TMP"/all_but_cortex.nii.gz -g 0 -n 1 -g 1 -n 0 -g 2 -n 1
AimsReplaceLevel -i "$TMP"/all_but_cortex.nii.gz -o "$TMP"/negative_outside_cortex.nii.gz -g 1 -n -1
AimsFileConvert -i "$TMP"/negative_outside_cortex.nii.gz -o "$TMP"/negative_outside_cortex_S32.nii.gz -t S32

AimsMerge -i "$TMP"/negative_outside_cortex_S32.nii.gz -M "$TMP"/CSF_labelled_interface.nii.gz -m sv -o "$TMP"/CSF_labelled_interface_negative_outside.nii.gz
AimsMerge -i "$TMP"/CSF_labelled_interface_negative_outside.nii.gz -M "$TMP"/white_labelled_interface.nii.gz -m ao -v 200000000 -o propvol_CSF_labels.nii.gz

AimsMerge -i "$TMP"/negative_outside_cortex_S32.nii.gz -M "$TMP"/white_labelled_interface.nii.gz -m sv -o "$TMP"/white_labelled_interface_negative_outside.nii.gz
AimsMerge -i "$TMP"/white_labelled_interface_negative_outside.nii.gz -M "$TMP"/CSF_labelled_interface.nii.gz -m ao -v 100000000 -o propvol_white_labels.nii.gz


ylPropagateAlongField --verbose \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --seeds ./propvol_CSF_labels.nii.gz \
    --step -0.01 \
    --target-label 200000000 \
    --output heat_CSF_labels_on_white.nii.gz
ylPropagateAlongField --verbose \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --seeds ./propvol_white_labels.nii.gz \
    --step 0.01 \
    --target-label 100000000 \
    --output heat_white_labels_on_CSF.nii.gz

python get_exchanged_propvol.py  # -> exchanged_propvol.nii.gz

enlarge_cortex_mask=false
if $enlarge_cortex_mask; then
    AimsMerge -i ./exchanged_propvol.nii.gz \
        -o ./exchanged_propvol_CSF_labels.nii.gz \
        -M ../classif_with_outer_boundaries.nii.gz \
        -m oo -l 15 -v 0
    AimsMerge -i ./exchanged_propvol.nii.gz \
        -o ./exchanged_propvol_white_labels.nii.gz \
        -M ../classif_with_outer_boundaries.nii.gz \
        -m oo -l 5 -v 0

    ylPropagateAlongField --verbose \
        --fieldx ../heat/heat_gradx.nii.gz \
        --fieldy ../heat/heat_grady.nii.gz \
        --fieldz ../heat/heat_gradz.nii.gz \
        --seeds ./exchanged_propvol_CSF_labels.nii.gz \
        --step -0.01 \
        --target-label 0 \
        --output "$TMP"/heat_CSF_on_bulk_raw.nii.gz \
        --dest-points ./heat_CSF_points_on_bulk.nii.gz
    ylPropagateAlongField --verbose \
        --fieldx ../heat/heat_gradx.nii.gz \
        --fieldy ../heat/heat_grady.nii.gz \
        --fieldz ../heat/heat_gradz.nii.gz \
        --seeds ./exchanged_propvol_white_labels.nii.gz \
        --step 0.01 \
        --target-label 0 \
        --output "$TMP"/heat_white_on_bulk_raw.nii.gz \
        --dest-points ./heat_white_points_on_bulk.nii.gz

    AimsReplaceLevel -i ../classif_with_outer_boundaries.nii.gz \
        -o "$TMP"/cortex_mask.nii.gz \
        -g 0 -n 0 -g 5 -n 1 -g 1 -n 1 -g 15 -n 1 -g 2 -n 0
else
    ylPropagateAlongField --verbose \
        --fieldx ../heat/heat_gradx.nii.gz \
        --fieldy ../heat/heat_grady.nii.gz \
        --fieldz ../heat/heat_gradz.nii.gz \
        --seeds ./exchanged_propvol.nii.gz \
        --step -0.01 \
        --target-label 0 \
        --output "$TMP"/heat_CSF_on_bulk_raw.nii.gz \
        --dest-points ./heat_CSF_points_on_bulk.nii.gz
    ylPropagateAlongField --verbose \
        --fieldx ../heat/heat_gradx.nii.gz \
        --fieldy ../heat/heat_grady.nii.gz \
        --fieldz ../heat/heat_gradz.nii.gz \
        --seeds ./exchanged_propvol.nii.gz \
        --step 0.01 \
        --target-label 0 \
        --output "$TMP"/heat_white_on_bulk_raw.nii.gz \
        --dest-points ./heat_white_points_on_bulk.nii.gz

    AimsReplaceLevel -i ../classif_with_outer_boundaries.nii.gz \
        -o "$TMP"/cortex_mask.nii.gz \
        -g 0 -n 0 -g 5 -n 0 -g 1 -n 1 -g 15 -n 0 -g 2 -n 0
fi

AimsMask -i "$TMP"/heat_CSF_on_bulk_raw.nii.gz \
    -o ./heat_CSF_on_bulk.nii.gz \
    -m "$TMP"/cortex_mask.nii.gz
AimsMask -i "$TMP"/heat_white_on_bulk_raw.nii.gz \
    -o ./heat_white_on_bulk.nii.gz \
    -m "$TMP"/cortex_mask.nii.gz

python relabel_conjunction.py  # -> ./conjunction.nii.gz

#AimsConnectComp -i conjunction.nii.gz -o conjunction_connected.nii.gz

ylMergeCortexColumnRegions \
    -i ./conjunction.nii.gz \
    -o merged.nii.gz \
    --proj-csf ./heat_CSF_points_on_bulk.nii.gz \
    --proj-white ./heat_white_points_on_bulk.nii.gz \
    --goal-diametre 0.5 \
    --max-thickness 2 \
    --verbose 2 \
&& python relabel.py \
&& python randomize_labels.py \
# && ylMakeRegionQualityMap \
#     -i ./merged.nii.gz \
#     -o region-roundness.nii.gz \
#     --proj-csf ./heat_CSF_points_on_bulk.nii.gz \
#     --proj-white ./heat_white_points_on_bulk.nii.gz \
#     --roundness-weight 1 \
#     --fullness-weight 0 \
#     --extension-weight 0 \
# && ylMakeRegionQualityMap \
#     -i ./merged.nii.gz \
#     -o region-fullness.nii.gz \
#     --proj-csf ./heat_CSF_points_on_bulk.nii.gz \
#     --proj-white ./heat_white_points_on_bulk.nii.gz \
#     --roundness-weight 0 \
#     --fullness-weight 1 \
#     --extension-weight 0 \
# && ylMakeRegionQualityMap \
#     -i ./merged.nii.gz \
#     -o region-extension.nii.gz \
#     --proj-csf ./heat_CSF_points_on_bulk.nii.gz \
#     --proj-white ./heat_white_points_on_bulk.nii.gz \
#     --roundness-weight 0 \
#     --fullness-weight 0 \
#     --extension-weight 1 \
&& AimsFileConvert \
    -i merged_randomized.nii.gz \
    -o merged_randomized_S16.nii.gz \
    -t S16 \
&& AimsGraphConvert \
    -i merged_randomized_S16.nii.gz \
    -o merged_randomized.arg \
    --bucket


N=1700 ; python \
    /volatile/yl232319/bv/src/perso/leprince/tools/proj_point_cloud.py \
    heat_CSF_points_on_bulk.nii.gz \
    merged_randomized.nii.gz \
    csf_projected_points.mesh \
    $N \
&& python \
    /volatile/yl232319/bv/src/perso/leprince/tools/proj_point_cloud.py \
    heat_white_points_on_bulk.nii.gz \
    merged_randomized.nii.gz \
    white_projected_points.mesh \
    $N \
&& AimsThreshold \
    -i merged_randomized.nii.gz \
    -b -m eq -t $N \
    -o focused_region.nii.gz \
&& AimsGraphConvert --bucket \
    -i focused_region.nii.gz \
    -o focused_region.arg


## Does not work with S32 labels (AimsVoronoi)
# AimsReplaceLevel -i heat_CSF_labels_on_white.nii.gz -o "$TMP"/heat_CSF_labels_on_white_negative_inside.nii.gz -g 0 -n -1
# AimsVoronoi -i "$TMP"/heat_CSF_labels_on_white_negative_inside.nii.gz -o "$TMP"/heat_CSF_labels_on_white_negative_inside_no_holes.nii.gz -d 200000000 -f -1
# AimsMerge -i "$TMP"/heat_CSF_labels_on_white_negative_inside_no_holes.nii.gz -M heat_CSF_labels_on_white.nii.gz -m oo -l 0 -v 0 -o heat_CSF_labels_on_white_no_holes.nii.gz


# Problem: gradient unreliable near boundaries due to huge values
ylPropagateAlongField --verbose \
    --fieldx ../dist/distmap_from_white_gradx.nii.gz \
    --fieldy ../dist/distmap_from_white_grady.nii.gz \
    --fieldz ../dist/distmap_from_white_gradz.nii.gz \
    --seeds ./propvol_CSF_labels.nii.gz \
    --step 0.01 \
    --target-label 200000000 \
    --output dist_CSF_labels_on_white.nii.gz
ylPropagateAlongField --verbose \
    --fieldx ../dist/distmap_from_CSF_gradx.nii.gz \
    --fieldy ../dist/distmap_from_CSF_grady.nii.gz \
    --fieldz ../dist/distmap_from_CSF_gradz.nii.gz \
    --seeds ./propvol_white_labels.nii.gz \
    --step 0.01 \
    --output dist_white_labels_on_CSF.nii.gz
