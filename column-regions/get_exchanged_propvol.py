# -*- coding: utf-8; -*-
import numpy as np
from soma import aims

CSF_labels_on_white = aims.read("heat_CSF_labels_on_white.nii.gz")
white_labels_on_CSF = aims.read("heat_white_labels_on_CSF.nii.gz")
classif = aims.read("../classif.nii.gz")
output = aims.Volume(CSF_labels_on_white)

np_CSF_labels_on_white = np.asarray(CSF_labels_on_white)
np_white_labels_on_CSF = np.asarray(white_labels_on_CSF)
np_classif = np.asarray(classif)
np_output = np.asarray(output)

white_mask = (np_classif == 2)
CSF_mask = (np_classif == 0)

np_output[white_mask] = np_CSF_labels_on_white[white_mask]
np_output[CSF_mask] = np_white_labels_on_CSF[CSF_mask]

aims.write(output, "raw_exchanged_labels.nii.gz")


# These “failed components” will probably be separated by connexity
#AimsReplaceLevel -i raw_exchanged_labels.nii.gz -o exchanged_labels.nii.gz -g 100000000 -n 0 -g 200000000 -n 0


import subprocess
subprocess.check_call(["AimsConnectComp",
                       "-i",
                       "raw_exchanged_labels.nii.gz",
                       "-o",
                       "connected_exchanged_labels.nii.gz"])


# The background is cut in one big region + many small, restore it then relabel
propvol = aims.read("connected_exchanged_labels.nii.gz")
np_propvol = np.asarray(propvol)
exclusion_mask = (np_CSF_labels_on_white == -1)
bulk_mask = (np_CSF_labels_on_white == 0)
np_propvol[bulk_mask] = 0
np_propvol[exclusion_mask] = -1

def relabel_positive_labels(volume):
    size_x = volume.getSizeX()
    size_y = volume.getSizeY()
    size_z = volume.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in xrange(size_z):
        for y in xrange(size_y):
            for x in xrange(size_x):
                old_label = volume.at(x, y, z)
                if old_label > 0:
                    try:
                        new_label = old_to_new_labels[old_label]
                    except KeyError:
                        new_label = next_label
                        old_to_new_labels[old_label] = new_label
                        next_label += 1
                    volume.setValue(new_label, x, y, z)

relabel_positive_labels(propvol)
aims.write(propvol, "exchanged_propvol.nii.gz")
