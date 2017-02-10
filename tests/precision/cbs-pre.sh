#! /bin/sh -e

mkdir CBS
mkdir CBS/Equidistant
mkdir CBS/Equivolumic
mkdir CBS/tmp

layout_template=$(dirname "$0")/evaluate.LayoutXML.in
sed -e "s~EVALUATION_DIR~${PWD}/CBS~g" \
    < "$layout_template" > CBS/evaluate.LayoutXML

AimsThreshold -b --fg 1 -m eq -t 200 \
    -i classif.nii.gz -o CBS/white_proba.nii.gz
AimsThreshold -b --fg 1 -m di -t 0 \
    -i classif.nii.gz -o CBS/not_CSF_proba.nii.gz
