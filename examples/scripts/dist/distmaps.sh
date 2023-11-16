#! /bin/sh -e
#
# Copyright Forschungszentrum Jülich GmbH (2018).
# Copyright CEA (2014).
# Copyright Université Paris XI (2014).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved. This file is offered as-is,
# without any warranty.

python -m capsul run highres_cortex.capsul.processes.Distmaps \
    classif=../classif.nii.gz \
    distwhite=distwhite.nii.gz \
    distCSF=distCSF.nii.gz \
    classif_with_outer_boundaries=../classif_with_outer_boundaries.nii.gz
