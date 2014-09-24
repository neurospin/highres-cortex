# coding: utf-8
from soma import aims
from soma import aimsalgo

import sys
n_iter = int(sys.argv[1])
time_step = float(sys.argv[2])

heatmap_before = aims.read("./heat.nii.gz", 1)

mask = aims.read("./all_but_cortex.nii", 1)
aimsmask = aims.AimsData(mask)  # Important for reference-counting!

diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step)
diff.setMask(aimsmask, 32767)
heatmap = diff.doSmoothing(heatmap_before, n_iter, True)

aims.write(heatmap, "./heat.nii.gz")
