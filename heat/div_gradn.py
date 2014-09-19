import numpy as np
from soma import aims

heatmap_volume = aims.read("heat.nii.gz")
header = heatmap_volume.header()
voxel_size_x, voxel_size_y, voxel_size_z = header["voxel_size"][:3]

heat = np.asarray(heatmap_volume)

class HalfcenteredGradients:
    def __init__(self, array):
        self.ppp = array[1:, 1:, 1:]
        self.ppm = array[1:, 1:, :-1]
        self.pmp = array[1:, :-1, 1:]
        self.pmm = array[1:, :-1, :-1]
        self.mpp = array[:-1, 1:, 1:]
        self.mpm = array[:-1, 1:, :-1]
        self.mmp = array[:-1, :-1, 1:]
        self.mmm = array[:-1, :-1, :-1]

    def gradx(self):
        return (0.25 / voxel_size_x) * (self.ppp + self.pmp + self.ppm + self.pmm
                                        - (self.mpp + self.mmp + self.mpm + self.mmm))
    def grady(self):
        return (0.25 / voxel_size_y) * (self.ppp + self.mpp + self.ppm + self.mpm
                                        - (self.pmp + self.mmp + self.pmm + self.mmm))
    def gradz(self):
        return (0.25 / voxel_size_z) * (self.ppp + self.mpp + self.pmp + self.mmp
                                        - (self.ppm + self.mpm + self.pmm + self.mmm))
    def gradxyz(self):
        return self.gradx(), self.grady(), self.gradz()

grad = HalfcenteredGradients(heat)
gradx, grady, gradz = grad.gradxyz()
gradn = np.sqrt(gradx ** 2 + grady ** 2 + gradz ** 2)
with np.errstate(divide="ignore", invalid="ignore"):
    gradx /= gradn
    grady /= gradn
    gradz /= gradn

grad = HalfcenteredGradients(gradx)
ggradx = grad.gradx()
grad = HalfcenteredGradients(grady)
ggrady = grad.grady()
grad = HalfcenteredGradients(gradz)
ggradz = grad.gradz()

div = ggradx + ggrady + ggradz

div_volume = aims.Volume(heatmap_volume)
div_volume.fill(float("NaN"))
div_array = np.asarray(div_volume)

div_array[1:-1, 1:-1, 1:-1] = div

aims.write(div_volume, "heat_div_gradn.nii.gz")
