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

from __future__ import absolute_import, division, print_function

import sys

from soma import aims

from six.moves import range


def relabel(labels):
    output = aims.Volume(labels)
    size_x = output.getSizeX()
    size_y = output.getSizeY()
    size_z = output.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in range(size_z):
        for y in range(size_y):
            for x in range(size_x):
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


def relabel_files(input_filename, output_filename):
    input_vol = aims.read(input_filename)
    output_vol = relabel(input_vol)
    aims.write(output_vol, output_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    import argparse
    parser = argparse.ArgumentParser(
        description="""\
Assign new consecutive labels to an existing label image
""")
    parser.add_argument("input", help="input label image")
    parser.add_argument("output", help="output label image")

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return relabel_files(
        args.input,
        args.output) or 0


if __name__ == "__main__":
    sys.exit(main())
