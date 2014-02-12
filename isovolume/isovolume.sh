AimsThreshold -b --fg 1 -m eq -t 100 \
    -i ../classif.nii.gz \
    -o ./domain.nii.gz

ylAdvectTubes --verbose \
    --step 0.05 \
    --domain ./domain.nii.gz \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --divergence ../heat/heat_div_gradn.nii.gz \
    --output-volumes white-tube-volumes.nii.gz \
    --output-surfaces white-tube-surfaces.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i white-tube-volumes.nii.gz \
    -i white-tube-surfaces.nii.gz \
    -o white-tube-VoverS.nii.gz

ylAdvectTubes --verbose \
    --step -0.05 \
    --domain ./domain.nii.gz \
    --fieldx ../heat/heat_gradx.nii.gz \
    --fieldy ../heat/heat_grady.nii.gz \
    --fieldz ../heat/heat_gradz.nii.gz \
    --divergence ../heat/heat_div_gradn.nii.gz \
    --output-volumes pial-tube-volumes.nii.gz \
    --output-surfaces pial-tube-surfaces.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i pial-tube-volumes.nii.gz \
    -i pial-tube-surfaces.nii.gz \
    -o pial-tube-VoverS.nii.gz

cartoLinearComb.py -f 'I1+I2' \
    -i pial-tube-volumes.nii.gz \
    -i white-tube-volumes.nii.gz \
    -o total-tube-volumes.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i pial-tube-volumes.nii.gz \
    -i total-tube-volumes.nii.gz \
    -o pial-volume-fraction.nii.gz
