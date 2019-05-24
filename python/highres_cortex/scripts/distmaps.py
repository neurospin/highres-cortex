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

import highres_cortex.cortex_topo


def compute_distmaps_files(classif_filename, output_distwhite_filename,
                           output_distCSF_filename, output_classif_filename):
    classif = aims.read(classif_filename)

    dist_from_white = highres_cortex.cortex_topo.signed_distance(
        classif, [100], [200], 150)
    aims.write(dist_from_white, output_distwhite_filename)

    dist_from_CSF = highres_cortex.cortex_topo.signed_distance(
        classif, [100], [0], 50)
    aims.write(dist_from_CSF, output_distCSF_filename)

    aims.write(classif, output_classif_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    import argparse
    parser = argparse.ArgumentParser(
        description="""\
Compute the signed distance to white matter and to CSF
""")
    parser.add_argument("classif", help="classification image of the cortex "
                        "(100 inside, 0 in CSF, 200 in white matter)")
    parser.add_argument("output_distwhite", help="signed Euclidean distance "
                        "to the white matter boundary")
    parser.add_argument("output_distCSF", help="signed Euclidean distance "
                        "to the CSF boundary")
    parser.add_argument("output_classif_with_boundaries",
                        help="classification image of the cortex (100 inside, "
                        "0 in CSF, 200 in white matter, 50 on the CSF "
                        "boundary, 150 on the white matter boundary)")

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return compute_distmaps_files(
        args.classif,
        args.output_distwhite,
        args.output_distCSF,
        args.output_classif_with_boundaries) or 0


if __name__ == "__main__":
    sys.exit(main())
