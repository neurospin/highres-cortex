import numpy as np
from soma import aims, aimsalgo

import yl.cortex_topo

classif = aims.read("../classif.nii")

dist_from_white = yl.cortex_topo.fastmarching_negative(
    classif, [100], [200], 150)
aims.write(dist_from_white, "./distwhite.nii.gz")

dist_from_CSF = yl.cortex_topo.fastmarching_negative(
    classif, [100], [0], 50)
aims.write(dist_from_CSF, "./distCSF.nii.gz")

aims.write(classif, "./classif_with_outer_boundaries.nii.gz")
