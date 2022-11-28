# -*- coding: utf-8 -*-
#
# Copyright CEA (2019).
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

import numpy as np
from soma import aims

from highres_cortex.cortex_topo import CSF_LABEL, WHITE_LABEL


def postprocess_equivolumetric_depth(input_filename, classif_filename,
                                     output_filename):
    depth_vol = aims.read(input_filename)
    classif_vol = aims.read(classif_filename)

    depth_arr = np.asarray(depth_vol)
    classif_arr = np.asarray(classif_vol)

    depth_arr[classif_arr == CSF_LABEL] = 0.0
    depth_arr[classif_arr == WHITE_LABEL] = 1.0

    header = depth_vol.header()
    header['cal_min'] = 0.0
    header['cal_max'] = 1.0
    header['intent_code'] = 1001  # NIFTI_INTENT_ESTIMATE
    header['intent_name'] = 'Equivol. depth'
    header['descrip'] = (
        'Equivolumetric cortical depth computed with highres-cortex'
    )

    aims.write(depth_vol, output_filename)


def parse_command_line(argv=sys.argv):
    """Parse the script's command line."""
    import argparse
    parser = argparse.ArgumentParser(
        description="""\
Post-process an equivolumetric depth image.

    - Set the outside of the brain (CSF) to 0.0
    - Set the white matter to 1.0
    - Set various Nifti header fields
""")
    parser.add_argument("input_image",
                        help="input image of equivolumetric depth")
    parser.add_argument("classif", help="classification image of the cortex "
                        "(100 inside, 0 in CSF, 200 in white matter)")
    parser.add_argument("output_image")

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    args = parse_command_line(argv)
    return postprocess_equivolumetric_depth(
        args.input_image,
        args.classif,
        args.output_image) or 0


if __name__ == "__main__":
    sys.exit(main())
