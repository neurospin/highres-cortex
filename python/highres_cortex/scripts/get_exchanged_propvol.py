#! /usr/bin/env python
# -*- coding: utf-8 -*-
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

import numpy as np
from soma import aims
import subprocess

def getExchangedPropagationVolume(CSF_labels_on_white, white_labels_on_CSF, classif, resultDir, keyWord):
    output = aims.Volume(CSF_labels_on_white)
    np_CSF_labels_on_white = np.asarray(CSF_labels_on_white)
    np_white_labels_on_CSF = np.asarray(white_labels_on_CSF)
    np_classif = np.asarray(classif)
    np_output = np.asarray(output)

    white_mask = (np_classif == 150)
    CSF_mask = (np_classif == 50)

    np_output[white_mask] = np_CSF_labels_on_white[white_mask]
    np_output[CSF_mask] = np_white_labels_on_CSF[CSF_mask]

    aims.write(output, resultDir + 'raw_exchanged_labels_%s.nii' %(keyWord))


    # These “failed components” will probably be separated by connexity
    #AimsReplaceLevel -i raw_exchanged_labels.nii.gz -o exchanged_labels.nii.gz -g 100000000 -n 0 -g 200000000 -n 0

    subprocess.check_call(["AimsConnectComp",
                        "-i", resultDir + 'raw_exchanged_labels_%s.nii' %(keyWord),
                        "-o", resultDir + 'connected_exchanged_labels_%s.nii' %(keyWord)])


    # The background is cut in one big region + many small, restore it then relabel
    propvol = aims.read(resultDir + 'connected_exchanged_labels_%s.nii' %(keyWord))
    np_propvol = np.asarray(propvol)
    exclusion_mask = (np_CSF_labels_on_white == -1)
    bulk_mask = (np_CSF_labels_on_white == 0)
    np_propvol[bulk_mask] = 0
    np_propvol[exclusion_mask] = -1

    def relabel_positive_labels(volume):
        size_x = volume.getSizeX()
        size_y = volume.getSizeY()
        size_z = volume.getSizeZ()
        old_to_new_labels = {}
        next_label = 1
        for z in xrange(size_z):
            for y in xrange(size_y):
                for x in xrange(size_x):
                    old_label = volume.at(x, y, z)
                    if old_label > 0:
                        try:
                            new_label = old_to_new_labels[old_label]
                        except KeyError:
                            new_label = next_label
                            old_to_new_labels[old_label] = new_label
                            next_label += 1
                        volume.setValue(new_label, x, y, z)

    relabel_positive_labels(propvol)
    return(propvol)


if __name__ == '__main__':
    heat_CSF = None
    heat_white = None
    classifFile = None
    resultDir = None
    keyWord = None

    parser = OptionParser('Get exchanged propagation volume -YL')
    parser.add_option('-s', dest='heat_CSF', help='heat_CSF_labels_on_white')   
    parser.add_option('-w', dest='heat_white', help='heat_white_labels_on_CSF') 
    parser.add_option('-c', dest='classifFile', help='classif_with_outer_boundaries')
    parser.add_option('-d', dest='resultDir', help='directory for results')
    parser.add_option('-k', dest='keyWord', help='keyword for results')

    options, args = parser.parse_args(sys.argv)
    print options
    print args

    if options.heat_CSF is None:
        print >> sys.stderr, 'New: exit. no heat_CSF_labels_on_white given'
        sys.exit(1)
    else:
        heat_CSF = options.heat_CSF
            
    if options.heat_white is None:
        print >> sys.stderr, 'New: exit. no heat_white_labels_on_CSF given'
        sys.exit(1)
    else:
        heat_white = options.heat_white
        
    if options.classifFile is None:
        print >> sys.stderr, 'New: exit. no classification file given'
        sys.exit(1)
    else:
        classifFile = options.classifFile      
 
    if options.resultDir is None:
        print >> sys.stderr, 'New: exit. no directory for results given'
        sys.exit(1)
    else:
        resultDir = options.resultDir    
        
    if options.keyWord is None:
        print >> sys.stderr, 'New: exit. no keyWord given'
        sys.exit(1)
    else:
        keyWord = options.keyWord      
     
      
    CSF_labels_on_white = aims.read(heat_CSF)
    white_labels_on_CSF = aims.read(heat_white)    
    classif = aims.read(classifFile)
    output = aims.Volume(CSF_labels_on_white)

    exchangedPropVol = getExchangedPropagationVolume(CSF_labels_on_white, white_labels_on_CSF, classif, resultDir, keyWord)
    aims.write(exchangedPropVol, resultDir + "exchanged_propvol_%s.nii.gz" %(keyWord))
    
    
    
