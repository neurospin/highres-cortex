# -*- coding: utf-8; -*-
import numpy as np
from soma import aims

input_labels = aims.read("./merged.nii.gz")

def relabel(labels):
    output = aims.Volume(labels)
    size_x = output.getSizeX()
    size_y = output.getSizeY()
    size_z = output.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in xrange(size_z):
        for y in xrange(size_y):
            for x in xrange(size_x):
                label = labels.at(x, y, z)
                if label == 0:
                    new_label = 0
                else:
                    try:
                        new_label = old_to_new_labels[label]
                    except KeyError:
                        new_label = next_label
                        old_to_new_labels[label] = new_label
                        next_label += 1
                output.setValue(new_label, x, y, z)
    return output

output = relabel(input_labels)
aims.write(output, "merged_relabelled.nii")
