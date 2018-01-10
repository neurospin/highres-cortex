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

import sys
import numpy as np
from soma import aims

def relabel_conjunctions(labels1, labels2):
    output = aims.Volume(labels1)
    output.fill(0)
    size_x = output.getSizeX()
    size_y = output.getSizeY()
    size_z = output.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in xrange(size_z):
        for y in xrange(size_y):
            for x in xrange(size_x):
                labels = (labels1.at(x, y, z), labels2.at(x, y, z))
                # Negative means outside propagation region
                if labels[0] < 0 or labels[1] < 0:
                    continue
                # Zeros are failed propagations, they should not be aggregated
                # together
                if labels[0] == 0 or labels[1] == 0:
                    new_label = next_label
                    next_label += 1
                else:
                    try:
                        new_label = old_to_new_labels[labels]
                    except KeyError:
                        new_label = next_label
                        old_to_new_labels[labels] = new_label
                        next_label += 1
                output.setValue(new_label, x, y, z)
    sys.stderr.write("{0}: {1} regions in conjunction\n"
                     .format(sys.argv[0], next_label - 1))
    return output

    
if __name__ == '__main__':
    CSF_labelsFile = None
    white_labelsFile = None
    resultDir = None
    keyWord = None

    parser = OptionParser('Get exchanged propagation volume -YL')
    parser.add_option('-s', dest='CSF_labelsFile', help='heat_CSF_labels_on_white')   
    parser.add_option('-w', dest='white_labelsFile', help='heat_white_labels_on_CSF') 
    parser.add_option('-d', dest='resultDir', help='directory for results')
    parser.add_option('-k', dest='keyWord', help='keyword for results')

    options, args = parser.parse_args(sys.argv)
    print options
    print args

    if options.CSF_labelsFile is None:
        print >> sys.stderr, 'New: exit. no heat_CSF_labels_on_white given'
        sys.exit(1)
    else:
        CSF_labelsFile = options.CSF_labelsFile
            
    if options.white_labelsFile is None:
        print >> sys.stderr, 'New: exit. no heat_white_labels_on_CSF given'
        sys.exit(1)
    else:
        white_labelsFile = options.white_labelsFile     
 
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
     

     
     
    CSF_labels = aims.read(CSF_labelsFile)
    white_labels = aims.read(white_labelsFile)
    output = relabel_conjunctions(CSF_labels, white_labels)
    aims.write(output, resultDir + 'conjunction_%s.nii.gz' %(keyWord))
    
    
    
    
    
    
    
    
