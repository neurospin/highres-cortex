#!/bin/sh -e

AimsThreshold -b -m eq -t 100 \
    -i ../classif.nii.gz \
    -o cortex_only.nii
AimsDilation -e 6 \
    -i cortex_only.nii \
    -o dilated_cortex.nii
AimsMerge -m oo -l 0 -v -1 \
    -i ../classif.nii.gz \
    -M dilated_cortex.nii \
    -o classif_near_cortex.nii

python distmaps.py

AimsThreshold -m di -t -1 \
    -i ./classif_with_outer_boundaries.nii \
    -o classif_with_background.nii
AimsMerge -m sv \
    -i ../classif.nii.gz \
    -M classif_with_background.nii \
    -o ../classif_with_outer_boundaries.nii.gz
