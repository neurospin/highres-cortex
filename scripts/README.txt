Original data: MR image (space.nii.gz) and its voxel-wise segmentation
(classif.nii.gz).

By convention, uninteresting intermediary files are in uncompressed NIfTI
format (.nii), while interesting files are in compressed NIfTI (.nii.gz).

Each script should be run within its containing directory (i.e. cd dist).

1. Run dist/distmaps.sh.

2. Run heat/heat.sh. Check full relaxation of the heat propagation.

3. Run isovolume/isovolume.sh. The isovolume depth metric is output in
   pial-volume-fraction.nii.gz.

4. Run column-regions/column-regions.sh. The output is in
   merged_randomized.nii.gz. The main parameter is --goal-diametre at the end
   of the script, it controls the typical diameter of merged regions (in
   millimetres).
