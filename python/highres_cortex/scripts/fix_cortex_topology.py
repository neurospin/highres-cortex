#! /usr/bin/env python

import argparse
import sys

from soma import aims

import highres_cortex.cortex_topo


def fix_cortex_topology_files(input_filename, output_filename,
                              filling_size, fclosing):
    """Call highres_cortex.cortex_topo.fix_cortex_topology on files."""
    input_volume = aims.read(input_filename)
    output = highres_cortex.cortex_topo.fix_cortex_topology(
        input_volume, filling_size, fclosing)
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
