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

import random
import sys

from soma import aims


def randomize_labels(labels):
    import numpy as np
    np_input_labels = np.asarray(labels)
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


def randomize_labels_files(input_filename, output_filename):
    input_vol = aims.read(input_filename)
    output_vol = randomize_labels(input_vol)
    aims.write(output_vol, output_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    import argparse
    parser = argparse.ArgumentParser(
        description="""\
Randomize the labels of an image with contiguous labels
""")
    parser.add_argument("input")
    parser.add_argument("output")

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return randomize_labels_files(
        args.input,
        args.output) or 0


if __name__ == "__main__":
    sys.exit(main())
