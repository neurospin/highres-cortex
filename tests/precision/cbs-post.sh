#! /bin/sh -e

cartoLinearComb.py -f '1-I1' \
    -i CBS/Equivolumic/_layering.nii.gz \
    -o CBS/Equivolumic/inverted_layering.nii.gz
cartoLinearComb.py -f '1-I1' \
    -i CBS/Equidistant/_layering.nii.gz \
    -o CBS/Equidistant/inverted_layering.nii.gz
