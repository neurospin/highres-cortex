#! /usr/bin/env python

from soma import aims
import numpy

inner_radius = 2
outer_radius = 5

voxel_size = 0.5

size_in_voxels = int(2 * (outer_radius / voxel_size) + 1)

classif = aims.Volume(size_in_voxels, size_in_voxels, size_in_voxels,
                      dtype="S16")
classif.header()["voxel_size"] = [voxel_size, voxel_size, voxel_size, 1]

classif.fill(0)

np_classif = numpy.asarray(classif)

s = slice(-voxel_size * (size_in_voxels // 2),
          voxel_size * (size_in_voxels // 2),
          size_in_voxels * 1j)
xgrid, ygrid, zgrid = numpy.ogrid[s, s, s]

np_classif[xgrid ** 2 + ygrid ** 2 + zgrid ** 2 < outer_radius ** 2] = 100
np_classif[xgrid ** 2 + ygrid ** 2 + zgrid ** 2 < inner_radius ** 2] = 200

aims.write(classif, "classif.nii.gz")
