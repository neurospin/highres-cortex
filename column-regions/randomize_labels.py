# -*- coding: utf-8; -*-
import random
from soma import aims

input_labels = aims.read("./merged_relabelled.nii.gz")

def relabel(labels):
    import numpy as np
    np_input_labels = np.asarray(input_labels)
    max_label = np.max(np_input_labels)
    nonzero_labels = list(range(1, max_label + 1))
    random.shuffle(nonzero_labels)
    new_labels = [0] + nonzero_labels

    output = aims.Volume(labels)
    size_x = output.getSizeX()
    size_y = output.getSizeY()
    size_z = output.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in xrange(size_z):
        for y in xrange(size_y):
            for x in xrange(size_x):
                old_label = labels.at(x, y, z)
                if old_label >= 0:
                    new_label = new_labels[old_label]
                else:
                    new_label = 0
                output.setValue(new_label, x, y, z)
    return output

output = relabel(input_labels)
aims.write(output, "merged_randomized.nii.gz")
