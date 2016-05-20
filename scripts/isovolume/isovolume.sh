#! /bin/sh -e
#
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

AimsThreshold -b --fg 1 -m eq -t 100 \
    -i ../classif.nii.gz \
    -o ./domain.nii.gz

ylAdvectTubes --verbose \
    --step 0.05 \
    --domain ./domain.nii.gz \
    --grad-field ../heat/heat.nii.gz \
    --divergence ../heat/heat_div_gradn.nii.gz \
    --output-volumes white-tube-volumes.nii.gz \
    --output-surfaces white-tube-surfaces.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i white-tube-volumes.nii.gz \
    -i white-tube-surfaces.nii.gz \
    -o white-tube-VoverS.nii.gz

ylAdvectTubes --verbose \
    --step -0.05 \
    --domain ./domain.nii.gz \
    --grad-field ../heat/heat.nii.gz \
    --divergence ../heat/heat_div_gradn.nii.gz \
    --output-volumes pial-tube-volumes.nii.gz \
    --output-surfaces pial-tube-surfaces.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i pial-tube-volumes.nii.gz \
    -i pial-tube-surfaces.nii.gz \
    -o pial-tube-VoverS.nii.gz

cartoLinearComb.py -f 'I1+I2' \
    -i pial-tube-volumes.nii.gz \
    -i white-tube-volumes.nii.gz \
    -o total-tube-volumes.nii.gz
cartoLinearComb.py -f 'I1/I2' \
    -i pial-tube-volumes.nii.gz \
    -i total-tube-volumes.nii.gz \
    -o pial-volume-fraction.nii.gz

AimsMerge -m oo -l 200 -v 1 \
    -i pial-volume-fraction.nii.gz \
    -M ../classif.nii.gz \
    -o equivolumic_depth.nii.gz
