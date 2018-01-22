#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright Télécom ParisTech (2015).
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

import argparse
import sys

from soma import aims

import highres_cortex.cortex_topo


def fix_cortex_topology_files(input_filename, output_filename,
                              filling_size, fclosing):
    """Call highres_cortex.cortex_topo.fix_cortex_topology on files."""
    input_volume = aims.read(input_filename)

    try:
        output = highres_cortex.cortex_topo.fix_cortex_topology(
            input_volume, filling_size, fclosing)
    except OSError as exc:
        print("error: the VipHomotopic command cannot be"
              " found or executed ({0}). Please make sure that"
              " Morphologist is properly installed and that the"
              " command is in your PATH.".format(exc.strerror))
        return 1
    except subprocess.CalledProcessError as exc:
        print("error: the VipHomotopic command returned an error code ({0})."
              " Please inspect its output above for more information."
              .format(exc.returncode))
        return 1

    # BUG: aims.write offers no error checking, so the program will exit
    # successfully even if writing fails
    aims.write(output, output_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    parser = argparse.ArgumentParser(
        description="""\
Impose the topology of a hollow sphere onto the cortex in a voxelwise
segmentation, which uses the following labels: 100 in the cortex itself, 0
outside (CSF), 200 inside (white matter). In the output, the cortex is defined
using 6-connectivity, each other compartment using 26-connectivity.
""")
    parser.add_argument("input",
                        help="3D volume containing the input segmentation")
    parser.add_argument("output",
                        help="output 3D volume")
    parser.add_argument("--filling-size", type=float, default=2.,
                        help="""\
The size, in millimetres, of the largest holes in either cortical boundary that
will be filled. This must be smaller than the thinnest cortex in the image. The
default value is 2 mm, which is appropriate for a human brain.""")
    parser.add_argument("--fclosing", type=float, default=10.,
                        help="""\
The radius of the morphological closing which is used by VipHomotopic in
Cortical surface mode to retrieve the brain's outer envelope. The default
value, 10 mm, is appropriate for a human brain.""")

    args = parser.parse_args(argv[1:])
    if not args.filling_size >= 0:
        parser.error("filling_size must be a non-negative number")
    if not args.fclosing >= 0:
        parser.error("fclosing must be a non-negative number")
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return fix_cortex_topology_files(
        args.input, args.output, args.filling_size, args.fclosing) or 0


if __name__ == "__main__":
    sys.exit(main())
