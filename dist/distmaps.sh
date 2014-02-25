#!/bin/sh -e

python distmaps.py

AimsThreshold -m di -t -1 \
    -i ./classif_with_outer_boundaries.nii \
    -o classif_with_background.nii
AimsMerge -m sv \
    -i ../classif.nii.gz \
    -M classif_with_background.nii \
    -o ../classif_with_outer_boundaries.nii.gz
