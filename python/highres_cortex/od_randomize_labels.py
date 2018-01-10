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

import random
from soma import aims

def relabel(labels):
    import numpy as np
    np_input_labels = np.asarray(labels)
    max_label = np.max(np_input_labels)
    nonzero_labels = list(range(1, max_label + 1))
    random.shuffle(nonzero_labels)
    new_labels = [0] + nonzero_labels

    output = aims.Volume(labels)
    size_x = output.getSizeX()
    size_y = output.getSizeY()
    size_z = output.getSizeZ()
    old_to_new_labels = {}
    next_label = 1
    for z in xrange(size_z):
        for y in xrange(size_y):
            for x in xrange(size_x):
                old_label = labels.at(x, y, z)
                if old_label >= 0:
                    new_label = new_labels[old_label]
                else:
                    new_label = 0
                output.setValue(new_label, x, y, z)
    return output

    
if __name__ == '__main__':
    
    mergedFile = None
    resultDir = None
    keyWord = None

    parser = OptionParser('Get the randomized relabeled volume -YL')
    parser.add_option('-m', dest='mergedFile', help='mergedFile')   
    parser.add_option('-d', dest='resultDir', help='directory for results')
    parser.add_option('-k', dest='keyWord', help='keyword for results')

    options, args = parser.parse_args(sys.argv)
    print options
    print args

    if options.mergedFile is None:
        print >> sys.stderr, 'New: exit. no mergedFile given'
        sys.exit(1)
    else:
        mergedFile = options.mergedFile
             
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
    
    
    input_labels = aims.read(mergedFile)        
    output = relabel(input_labels)
    aims.write(output, resultDir + "merged_randomized_%s.nii.gz" %(keyWord))
