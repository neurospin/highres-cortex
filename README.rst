================
 highres-cortex
================

This is a collection of software designed to process 3D images of the cerebral cortex at a sub-millimetre scale, for example high-resolution MRI. In particular, it implements Bok’s equivolumetric model for extracting cortical layers while compensating for cortical curvature. If you use this work in an academic publication, **please cite** the relevant references (see also `<doc/references.bib>`_):

- Yann Leprince, Fabrice Poupon, Thierry Delzescaux, Dominique Hasboun, Cyril Poupon, et al.. *Combined Laplacian-equivolumic model for studying cortical lamination with ultra high field MRI (7 T)*. 2015 IEEE 12th International Symposium on Biomedical Imaging (ISBI), IEEE, Apr 2015, New York, United States. pp.580-583, DOI: `10.1109/ISBI.2015.7163940 <https://dx.doi.org/10.1109/ISBI.2015.7163940>`_.  https://hal-cea.archives-ouvertes.fr/cea-01119475

- Yann Leprince, Clara Fischer, Jean-François Mangin, Benoît Larrat, Sébastien Mériaux, et al.. *Architectonics-informed partition of the cortex at sub-millimetre resolution*. 20th Annual Meeting of the Organization for Human Brain Mapping (OHBM), Jun 2014, Hamburg, Germany. 5, pp.951, 2014, F1000Posters. https://hal-cea.archives-ouvertes.fr/cea-01074735



Installation
============

highres-cortex is released as part of the official BrainVISA containers since March 2021 (BrainVISA 5.0.0). Please refer to https://brainvisa.info/web/download.html for installation instructions.


Basic usage
===========

This package can be used on the command line, here is a short introduction. It is assumed that you are running the Singularity version on a Linux computer, it should work in the same way in the virtual machine (you just need to remove the `bv` prefix from all commands).

1. Prepare your input data: the input that is common to to all processes is ``classif``: a voxel-wise tissue classification image in signed 16-bit pixel type, with 0 for exterior voxels (CSF), 100 for cortical gray matter, and 200 for subcortical white matter.

2. Run the process that you are interested in. The common interface to all processes is the Capsul command-line, which you can call with ``bv python -m capsul``. Use ``bv python -m capsul --process-help <process_name>`` to get help for a specific process. Use ``bv python -m capsul <process_name> [parameter=value ...]`` to run a process.

   The most important processes are described below:

   - Equivolumetric depth according to Bok’s model can be computed with ``highres_cortex.capsul.isovolume``. The only mandatory input is ``classif``, the output is ``equivolumetric_depth``. You can fine-tune the process with optional parameters, most importantly ``advection_step_size`` can be adapted to the spatial resolution and required accuracy. For example::

         bv python -m capsul highres_cortex.capsul.isovolume classif=classif.nii.gz advection_step_size=0.03 equivolumetric_depth=equivolumetric_depth.nii.gz

   - Cortical thickness, according to the Laplace model, can be calculated with two different methods:

     - The upwinding method is very fast, and already has sub-pixel accurracy: ``highres_cortex.capsul.thickness_upw``. The only mandatory input is ``classif``, the output is ``thickness_image``.

     - The advection method is slower, but ``advection_step_size`` can be tuned for greater accuracy: ``highres_cortex.capsul.thickness_adv``.

   - For parcellating the cortex into volumetric traverses, ``highres_cortex.capsul.traverses`` can be used. The only mandatory input is ``classif``, the output is ``cortical_traverses``. The ``goal_diameter`` parameter controls the target diameter of merged regions (in millimetres). The ``advection_step_size`` parameter is also relevant for this process.

If you have used highres-cortex before the Capsul interface was introduced (beginning of 2018), you may be using the old shell scripts. See `<examples/scripts/>`_ for equivalent scripts that make use of the Capsul processes.


Contributing
============

This repository uses `pre-commit`_ to ensure that all committed code follows minimal quality standards. Please install it and configure it to run as a pre-commit hook in your local repository (note that this is done automatically by ``bv_maker``):

.. code-block:: shell

  # Install pre-commit in a virtual environment
  python3 -m venv venv/
  . venv/bin/activate
  pip install pre-commit

  pre-commit install  # install the pre-commit hook


Licence
=======

The source code of this work is placed under the CeCILL licence (see `<LICENCE.CeCILL.txt>`_). This library contains code that is under the GNU LGPL licence (see `<src/library/cortex_column_region_quality.tcc>`_), as a result, compiled code must be redistributed under the GNU General Public Licence (see `<LICENCE.GPLv3.txt>`_).

External code used in this repository
-------------------------------------

- Code for numerical diagonalization of 3×3 matrices (`<src/library/cortex_column_region_quality.tcc>`_) is Copyright 2006 Joachim Kopp, under the GNU LGPL v2.1 or later. Reference: Kopp, Joachim. ‘Efficient Numerical Diagonalization of Hermitian 3x3 Matrices’. *International Journal of Modern Physics C* 19, no. 03 (March 2008): 523–48. `arXiv:physics/0610206 <http://arxiv.org/abs/physics/0610206>`_.


.. _BrainVISA: http://brainvisa.info/
.. _Capsul: http://neurospin.github.io/capsul/
.. _Python: https://www.python.org/
.. _pre-commit: https://pre-commit.com/
.. _Introduction to bv_maker: https://bioproj.extra.cea.fr/redmine/projects/brainvisa-devel/wiki/How_to_compile_BrainVISA_projects
.. _BrainVISA download page: http://brainvisa.info/web/download.html

.. Copyright CEA (2014, 2015, 2021, 2022).
   Copyright Forschungszentrum Jülich GmbH (2016, 2017, 2018).
   Copyright Télécom ParisTech (2015, 2016).
   Copyright Université Paris XI (2014).

   Author: Yann Leprince <yann.leprince@cea.fr>.

   Copying and distribution of this file, with or without modification, are permitted in any medium without royalty provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.
