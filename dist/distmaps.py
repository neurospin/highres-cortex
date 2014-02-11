import numpy as np
from soma import aims, aimsalgo

import yl.cortex_topo

classif = aims.read("classif_near_cortex.nii")

dist_from_white = yl.cortex_topo.fastmarching_negative(
    classif, [100], [200], 150)
aims.write(dist_from_white, "./distwhite.nii.gz")

dist_from_CSF = yl.cortex_topo.fastmarching_negative(
    classif, [100], [0], 50)
aims.write(dist_from_CSF, "./distCSF.nii.gz")

aims.write(classif, "./classif_with_outer_boundaries.nii.gz")

gradient = aimsalgo.AimsGradient_FLOAT()
dist_from_white = aims.read("./distwhite.nii.gz")
gradX_from_white = gradient.X(dist_from_white)
gradY_from_white = gradient.Y(dist_from_white)
gradZ_from_white = gradient.Z(dist_from_white)
dist_from_CSF = aims.read("./distCSF.nii.gz")
gradX_from_CSF = gradient.X(dist_from_CSF)
gradY_from_CSF = gradient.Y(dist_from_CSF)
gradZ_from_CSF = gradient.Z(dist_from_CSF)
aims.write(gradX_from_white, "./distwhite_gradx.nii.gz")
aims.write(gradY_from_white, "./distwhite_grady.nii.gz")
aims.write(gradZ_from_white, "./distwhite_gradz.nii.gz")
aims.write(gradX_from_CSF, "./distCSF_gradx.nii.gz")
aims.write(gradY_from_CSF, "./distCSF_grady.nii.gz")
aims.write(gradZ_from_CSF, "./distCSF_gradz.nii.gz")

laplacian_from_white = gradient.X(gradX_from_white) + gradient.Y(gradY_from_white) + gradient.Z(gradZ_from_white)
laplacian_from_CSF = gradient.X(gradX_from_CSF) + gradient.Y(gradY_from_CSF) + gradient.Z(gradZ_from_CSF)
aims.write(laplacian_from_white, "./distwhite_laplacian.nii.gz")
aims.write(laplacian_from_CSF, "./distCSF_laplacian.nii.gz")
