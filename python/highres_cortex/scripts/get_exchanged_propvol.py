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

import os.path
import shutil
import sys
import tempfile
import six

import numpy as np
import soma.subprocess as subprocess
from soma import aims


def relabel_positive_labels(volume):
    size_x = volume.getSizeX()
    size_y = volume.getSizeY()
    size_z = volume.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in six.moves.xrange(size_z):
        for y in six.moves.xrange(size_y):
            for x in six.moves.xrange(size_x):
                old_label = volume.at(x, y, z)
                if old_label > 0:
                    try:
                        new_label = old_to_new_labels[old_label]
                    except KeyError:
                        new_label = next_label
                        old_to_new_labels[old_label] = new_label
                        next_label += 1
                    volume.setValue(new_label, x, y, z)


def get_exchanged_propvol_files(classif_filename,
                                CSF_labels_on_white_filename,
                                white_labels_on_CSF_filename,
                                output_filename):
    classif = aims.read(classif_filename)
    CSF_labels_on_white = aims.read(CSF_labels_on_white_filename)
    white_labels_on_CSF = aims.read(white_labels_on_CSF_filename)
    output = aims.Volume(CSF_labels_on_white)

    np_CSF_labels_on_white = np.asarray(CSF_labels_on_white)
    np_white_labels_on_CSF = np.asarray(white_labels_on_CSF)
    np_classif = np.asarray(classif)
    np_output = np.asarray(output)

    white_mask = (np_classif == 150)
    CSF_mask = (np_classif == 50)

    np_output[white_mask] = np_CSF_labels_on_white[white_mask]
    np_output[CSF_mask] = np_white_labels_on_CSF[CSF_mask]

    temp_dir = None
    try:
        temp_dir = tempfile.mkdtemp(prefix="hcortex")
        temp_filename = os.path.join(temp_dir, 'raw_exchanged_labels.nii')
        aims.write(output, temp_filename)

        # These “failed components” will probably be separated by connexity
        # AimsReplaceLevel -i raw_exchanged_labels.nii.gz \
        # -o exchanged_labels.nii.gz \
        # -g 100000000 -n 0 -g 200000000 -n 0

        subprocess.check_call(["AimsConnectComp",
                               "-i", "raw_exchanged_labels.nii",
                               "-o", "connected_exchanged_labels.nii"],
                              cwd=temp_dir)

        # The background is cut in one big region + many small, restore it then
        # relabel
        propvol = aims.read(
            os.path.join(temp_dir, "connected_exchanged_labels.nii"))
    finally:
        if temp_dir:
            shutil.rmtree(temp_dir)
    np_propvol = np.asarray(propvol)
    exclusion_mask = (np_CSF_labels_on_white == -1)
    bulk_mask = (np_CSF_labels_on_white == 0)
    np_propvol[bulk_mask] = 0
    np_propvol[exclusion_mask] = -1

    relabel_positive_labels(propvol)
    aims.write(propvol, output_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    import argparse
    parser = argparse.ArgumentParser(
        description="""\
Get exchanged propagation volume
""")
    parser.add_argument("classif_with_outer_boundaries", help="classification "
                        "image of the cortex (100 inside, 0 in CSF, 200 in "
                        "white matter, 50 on the CSF border, 150 on the white "
                        "matter border)")
    parser.add_argument("CSF_labels_on_white", help="labels of the CSF "
                        "projected onto the white matter boundary")
    parser.add_argument("white_labels_on_CSF", help="labels of the white "
                        "matter projected onto the CSF boundary")
    parser.add_argument("output", help="volume where each interface is "
                        "labelled with connected components facing the same "
                        "voxels of the other interface")

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return get_exchanged_propvol_files(
        args.classif_with_outer_boundaries,
        args.CSF_labels_on_white,
        args.white_labels_on_CSF,
        args.output) or 0


if __name__ == "__main__":
    sys.exit(main())
