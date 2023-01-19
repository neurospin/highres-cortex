#! /bin/sh -e
#
# Copyright Forschungszentrum Jülich GmbH (2018).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved. This file is offered as-is,
# without any warranty.

python -m capsul highres_cortex.capsul.isovolume \
    classif=../classif.nii.gz \
    verbosity=1 \
    advection_step_size=0.05 \
    equivolumetric_depth=equivolumic_depth.nii.gz
