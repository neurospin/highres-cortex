# coding: utf-8
from soma import aims
from soma import aimsalgo
classif = aims.read("../classif.nii.gz", 1)
conv = aims.Converter_Volume_S16_Volume_FLOAT()
classif_float = conv(classif)

mask = aims.read("./all_but_cortex.nii", 1)
aimsmask = aims.AimsData(mask)  # Important for reference-counting!

diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(0.1)
diff.setMask(aimsmask, 32767)
heatmap = diff.doSmoothing(classif_float, 500, True)

aims.write(heatmap, "./heat.nii.gz")


gradient = aimsalgo.AimsGradient_FLOAT()
gradx = gradient.X(heatmap)
grady = gradient.Y(heatmap)
gradz = gradient.Z(heatmap)
aims.write(gradx, "./heat_gradx.nii.gz")
aims.write(grady, "./heat_grady.nii.gz")
aims.write(gradz, "./heat_gradz.nii.gz")

laplacian = (gradient.X(gradx) +
             gradient.Y(grady) +
             gradient.Z(gradz))
aims.write(laplacian, "./heat_laplacian.nii.gz")
