#!/bin/sh -e

AimsThreshold -b -m di -t 100 \
    -i ../classif.nii.gz \
    -o ./all_but_cortex.nii
AimsFileConvert -t FLOAT \
    -i ../classif.nii.gz \
    -o heat.nii.gz

# Each run refines the previous one
python heat.py 500 0.01
python heat.py 100 0.001

# Normalized gradient's divergence
python div_gradn.py
