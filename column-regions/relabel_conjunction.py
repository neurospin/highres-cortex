# -*- coding: utf-8; -*-
import sys
import numpy as np
from soma import aims

CSF_labels = aims.read("./heat_CSF_on_bulk.nii.gz")
white_labels = aims.read("./heat_white_on_bulk.nii.gz")

def relabel_conjunctions(labels1, labels2):
    output = aims.Volume(labels1)
    output.fill(0)
    size_x = output.getSizeX()
    size_y = output.getSizeY()
    size_z = output.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in xrange(size_z):
        for y in xrange(size_y):
            for x in xrange(size_x):
                labels = (labels1.at(x, y, z), labels2.at(x, y, z))
                # Negative means outside propagation region
                if labels[0] < 0 or labels[1] < 0:
                    continue
                # Zeros are failed propagations, they should not be aggregated
                # together
                if labels[0] == 0 or labels[1] == 0:
                    new_label = next_label
                    next_label += 1
                else:
                    try:
                        new_label = old_to_new_labels[labels]
                    except KeyError:
                        new_label = next_label
                        old_to_new_labels[labels] = new_label
                        next_label += 1
                output.setValue(new_label, x, y, z)
    sys.stderr.write("{0}: {1} regions in conjunction\n"
                     .format(sys.argv[0], next_label - 1))
    return output

output = relabel_conjunctions(CSF_labels, white_labels)
aims.write(output, "conjunction.nii.gz")
