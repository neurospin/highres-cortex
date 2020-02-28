================
 highres-cortex
================

This is a collection of software designed to process 3D images of the cerebral cortex at a sub-millimetre scale, for example high-resolution MRI. In particular, it implements Bok’s equivolumetric model for extracting cortical layers while compensating for cortical curvature.

If you use this work in an academic publication, **please cite** the relevant references (see `<doc/references.bib>`_).


Basic usage
===========

This package can be used on the command line, here is a short introduction. See below for installation instructions.

1. Set up the necessary environment using ``bv_env.sh``::

     . </path/to/installation>/bin/bv_env.sh

2. Prepare your input data: the input that is common to to all processes is ``classif``: a voxel-wise tissue classification image in signed 16-bit pixel type, with 0 for exterior voxels (CSF), 100 for cortical gray matter, and 200 for subcortical white matter.

3. Run the process that you are interested in. The common interface to all processes is the Capsul command-line, which you can call with ``python -m capsul`` (or ``python -m capsul.run`` for Capsul >= 2.1.3). Use ``python -m capsul --process-help <process_name>`` to get help for a specific process. Use ``python -m capsul <process_name> [parameter=value ...]`` to run a process.

   The most important processes are described below:

   - Equivolumetric depth according to Bok’s model can be computed with ``highres_cortex.capsul.isovolume``. The only mandatory input is ``classif``, the output is ``equivolumetric_depth``. You can fine-tune the process with optional parameters, most importantly ``advection_step_size`` can be adapted to the spatial resolution and required accuracy. For example::

         python -m capsul highres_cortex.capsul.isovolume classif=classif.nii.gz advection_step_size=0.03 equivolumetric_depth=equivolumetric_depth.nii.gz

   - Cortical thickness, according to the Laplace model, can be calculated with two different methods:

     - The upwinding method is very fast, and already has sub-pixel accurracy: ``highres_cortex.capsul.thickness_upw``. The only mandatory input is ``classif``, the output is ``thickness_image``.

     - The advection method is slower, but ``advection_step_size`` can be tuned for greater accuracy: ``highres_cortex.capsul.thickness_adv``.

   - For parcellating the cortex into volumetric traverses, ``highres_cortex.capsul.traverses`` can be used. The only mandatory input is ``classif``, the output is ``cortical_traverses``. The ``goal_diameter`` parameter controls the target diameter of merged regions (in millimetres). The ``advection_step_size`` parameter is also relevant for this process.

If you have used highres-cortex before the Capsul interface was introduced (beginning of 2018), you may be using the old shell scripts. See `<examples/scripts/>`_ for equivalent scripts that make use of the Capsul processes.


Installation
============

Binary packages
---------------

We are planning on packaging highres-cortex, and making it available as part of the graphical installer of BrainVISA. For now you have to compile it as described below.


Automated compilation
---------------------

The script ``bootstrap_compile.sh`` can be used to download, configure, and build highres-cortex::

    wget https://github.com/neurospin/highres-cortex/raw/master/bootstrap_compile.sh
    chmod +x bootstrap_compile.sh
    ./bootstrap_compile.sh

The script is **interactive**, and will require some input from you.

- If you are on Ubuntu 14.04 or 16.04 LTS, this script will ensure that all the required dependencies are present on your machine, and propose their installation otherwise. For any other distribution, or any other version of Ubuntu, you have to install the dependencies manually.

- You will also be asked for a base directory to contain the downloaded sources, and the build process. You will need about **1.5 GB** of free space in this directory.

The step-by-step approach in the next section performs essentially the same steps as the bootstrapping script, use that if you want to customize the build process.


Step-by-step compilation
------------------------

You can compile this package as part of the BrainVISA_ source tree, which is based on CMake_ and uses a custom-made driver called ``bv_maker``.

1. Install the dependencies. Under Ubuntu, the required packages are: ``subversion git cmake make gcc g++ gfortran pkg-config libblitz0-dev libsigc++-2.0-dev libxml2-dev libqt4-dev libboost-dev zlib1g-dev libtiff-dev python2.7-dev python-sip-dev python-numpy python-six libqt4-sql-sqlite``.

2. Bootstrap the ``bv_maker`` tool::

     svn export --username brainvisa --password Soma2009 https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/branches/bug_fix /tmp/brainvisa-cmake
     cd /tmp/brainvisa-cmake
     cmake -DCMAKE_INSTALL_PREFIX=. .
     make install

   Detailed instructions can be found in this `Introduction to bv_maker`_ (login *brainvisa*, password *Soma2009*).

3. Create the configuration file for ``bv_maker`` at ``$HOME/.brainvisa/bv_maker.cfg``. Here is a minimal version of this file::

     [ source $HOME/brainvisa/source ]
       brainvisa brainvisa-cmake bug_fix
       brainvisa soma-base bug_fix
       brainvisa soma-io bug_fix
       brainvisa aims-free bug_fix
       brainvisa soma-workflow $CASA_BRANCH
       brainvisa capsul $CASA_BRANCH
       git https://github.com/neurospin/highres-cortex.git master highres-cortex

     [ build $HOME/brainvisa/build ]
       build_type = Release
       brainvisa-cmake bug_fix $HOME/brainvisa/source
       brainvisa-share bug_fix $HOME/brainvisa/source
       soma-base bug_fix $HOME/brainvisa/source
       soma-io bug_fix $HOME/brainvisa/source
       aims-free bug_fix $HOME/brainvisa/source
       + $HOME/brainvisa/source/highres-cortex

   Keep the following in mind if you want to customize this configuration file:
    - you need this line in the ``source`` section::

        git https://github.com/neurospin/highres-cortex.git master highres-cortex

    - you need this line in the ``build`` section::

        + </path/to/brainvisa/source>/highres-cortex

    - you need to enable the ``aims-free`` component and its dependencies ``brainvisa-cmake``, ``soma-base``, and ``soma-io``; alternatively, just enable the ``anatomist`` group, which is a superset of these.

4. Run ``/tmp/brainvisa-cmake/bin/bv_maker``, which will check out a local copy of the sources, configure them with cmake, and build thim with ``make``.

5. You can then run the software directly from ``$HOME/brainvisa/build``, as indicated in the `Basic usage`_ section.


Dependencies
============

- AIMS version 4.5 or later, an image processing library distributed as part of BrainVISA_.
- Boost_ version 1.49 or later.
- Python_ version 2.6 or later.
- CMake_ version 2.6 or later, with its extension ``brainvisa-cmake`` (distributed with BrainVISA_).
- Recommended: Capsul_ version 2 or later, used to combine the low-level building blocks into useful processing pipelines.
- Optional: the ``VipHomotopic`` command-line tool from the Morphologist image segmentation pipeline, distributed as a binary only tool with the BrainVISA_ installer.


Contributing
============

This repository uses `pre-commit`_ to ensure that all committed code follows minimal quality standards. Please install it and configure it to run as a pre-commit hook in your local repository:

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
.. _Boost: http://www.boost.org/
.. _CMake: http://www.cmake.org/
.. _Capsul: http://neurospin.github.io/capsul/
.. _Python: https://www.python.org/
.. _pre-commit: https://pre-commit.com/
.. _Introduction to bv_maker: https://bioproj.extra.cea.fr/redmine/projects/brainvisa-devel/wiki/How_to_compile_BrainVISA_projects
.. _BrainVISA download page: http://brainvisa.info/web/download.html

.. Copyright Forschungszentrum Jülich GmbH (2016, 2017, 2018).
   Copyright Télécom ParisTech (2015, 2016).
   Copyright CEA (2014, 2015).
   Copyright Université Paris XI (2014).

   Author: Yann Leprince <yann.leprince@ylep.fr>.

   Copying and distribution of this file, with or without modification, are permitted in any medium without royalty provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.
