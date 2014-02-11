#!/bin/sh -e

AimsThreshold -b -m di -t 100 \
    -i ../classif.nii.gz \
    -o ./all_but_cortex.nii.gz

python heat.py

# Normalized gradient
ylNumpyComb.py -f 'sqrt(I1**2 + I2**2 + I3**2)' \
    -i heat_gradx.nii.gz \
    -i heat_grady.nii.gz \
    -i heat_gradz.nii.gz \
    -o heat_grad_norm.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i heat_gradx.nii.gz \
    -i heat_grad_norm.nii.gz \
    -o heat_gradnx.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i heat_grady.nii.gz \
    -i heat_grad_norm.nii.gz \
    -o heat_gradny.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i heat_gradz.nii.gz \
    -i heat_grad_norm.nii.gz \
    -o heat_gradnz.nii.gz

python div_gradn.py

cartoLinearComb.py -f '-I1' \
    -i heat_gradx.nii.gz \
    -o heat_minus_gradx.nii.gz
cartoLinearComb.py -f '-I1' \
    -i heat_grady.nii.gz \
    -o heat_minus_grady.nii.gz
cartoLinearComb.py -f '-I1' \
    -i heat_gradz.nii.gz \
    -o heat_minus_gradz.nii.gz
cartoLinearComb.py -f '-I1' \
    -i heat_div_gradn.nii.gz \
    -o heat_minus_div_gradn.nii.gz
