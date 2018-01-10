#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright Forschungszentrum Jülich GmbH (2018).
# Copyright CEA (2014).
# Copyright Université Paris XI (2014).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# This file is part of highres-cortex, a collection of software designed
# to process high-resolution magnetic resonance images of the cerebral
# cortex.
#
# This software is governed by the CeCILL licence under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# licence as circulated by CEA, CNRS and INRIA at the following URL:
# <http://www.cecill.info/>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the licence, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of scientific
# software, that may mean that it is complicated to manipulate, and that
# also therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL licence and that you accept its terms.

import sys

from soma import aims


def relabel_conjunction(labels1, labels2):
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


def relabel_conjunction_files(labels1_filename, labels2_filename,
                              output_filename):
    labels1_vol = aims.read(labels1_filename)
    labels2_vol = aims.read(labels2_filename)
    output_vol = relabel_conjunction(labels1_vol, labels2_vol)
    aims.write(output_vol, output_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    import argparse
    parser = argparse.ArgumentParser(
        description="""\
Assign new labels to voxels that have the same pair of labels in both input images.
""")
    parser.add_argument("labels1")
    parser.add_argument("labels2")
    parser.add_argument("output")

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return relabel_conjunction_files(
        args.labels1,
        args.labels2,
        args.output) or 0

if __name__ == "__main__":
    sys.exit(main())
