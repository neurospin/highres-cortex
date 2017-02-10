================
 highres-cortex
================

This is a collection of software designed to process 3D images of the cerebral cortex at a sub-millimetre scale, for example high-resolution MRI. In particular, it implements Bok’s equivolumetric model for extracting cortical layers while compensating for cortical curvature.

If you use this work in an academic publication, **please cite** the relevant references (see `<doc/references.bib>`_).


Basic usage
===========

This package can be used on the command line, the main interface is through shell scripts. Here is a short introduction. See below for installation instructions.

1. Set up the necessary environment using ``bv_env.sh``::

     . </path/to/installation>/bin/bv_env.sh </path/to/installation>

2. Create a directory to process your data, and copy the contents of the `<scripts/>`_ sub-directory thereunder::

     cp -R </path/to>/scripts/* .

3. Put your input data into place. The only input is ``classif.nii.gz``, which should contain a voxel-wise tissue classification in signed 16-bit pixel type, with 0 for exterior voxels, 100 for cortical gray matter, and 200 for subcortical white matter.

4. Run the scripts that correspond to the output that you need. Each script should be run from within its own sub-directory (first ``cd`` to this directory, then run ``./script.sh``).

   - For the Laplace model, you need to run `heat.sh <scripts/heat/heat.sh>`_. This should work out of the box, without the need for tuning any parameter. Outputs are the Laplace field under ``heat/heat.nii.gz``, and the curvature field under ``heat_div_gradn.nii.gz``.

   - For Bok’s equivolumetric depth, you need to run `heat.sh <scripts/heat/heat.sh>`_, then `isovolume.sh <scripts/isovolume/isovolume.sh>`_. The result is output in ``isovolume/pial-volume-fraction.nii.gz``. You can tune parameters in ``isovolume.sh``, most importantly the step size (``--step``) can be adapted to the spatial resolution and required accuracy.

   - For calculating the Euclidian depth along Laplace traverses, you first need to run `heat.sh <scripts/heat/heat.sh>`_. Then, you have two choices:

     - Using the fast upwinding method in `upwind-euclidean.sh <scripts/upwind-euclidean/upwind-euclidean.sh>`_.

     - Using the slower advection method in `laplace-euclidean.sh <scripts/laplace-euclidean/laplace-euclidean.sh>`_.

   - For parcellating the cortex into volumetric traverses, you need to run `distmaps.sh <scripts/dist/distmaps.sh>`_, then `heat.sh <scripts/heat/heat.sh>`_, and finally `column-regions.sh <scripts/column-regions/column-regions.sh>`_. The main parameter is ``--goal-diameter`` at the end of the script, it controls the target diameter of merged regions (in millimetres). The ``--step`` parameter is also relevant here.


Installation
============

Binary packages
---------------

We are planning on packaging highres-cortex, and making it available as part of the graphical installer of BrainVISA. For now you have to compile it as described below.


Automated compilation
---------------------

The script ``bootstrap_compile.sh`` can be used to download, configure, and build highres-cortex::

    wget https://github.com/ylep/highres-cortex/raw/master/bootstrap_compile.sh
    chmod +x bootstrap_compile.sh
    ./bootstrap_compile.sh

The script is **interactive**, and will require some input from you.

- If you are on Ubuntu 14.04 or 16.04 LTS, this script will ensure that all the required dependencies are present on your machine, and propose their installation otherwise. For any other distribution, or any other version of Ubuntu, you have to install the dependencies manually.

- You will also be asked for a base directory to contain the downloaded sources, and the build process. You will need about **1.5 GB** of free space in this directory.

The step-by-step approach in the next section performs essentially the same steps as the bootstrapping script, use that if you want to customize the build process.


Step-by-step compilation
------------------------

You can compile this package as part of the BrainVISA_ source tree, which is based on CMake_ and uses a custom-made driver called ``bv_maker``.

1. Install the dependencies. Under Ubuntu, the required packages are: ``subversion git cmake make gcc g++ gfortran pkg-config libblitz0-dev libsigc++-2.0-dev libxml2-dev libqt4-dev libboost-dev zlib1g-dev libtiff-dev libgsl0-dev python2.7-dev python-sip-dev python-numpy python-six libqt4-sql-sqlite``.

2. Bootstrap the ``bv_maker`` tool::

     svn export --username brainvisa --password Soma2009 https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/branches/bug_fix /tmp/brainvisa-cmake
     cd /tmp/brainvisa-cmake
     cmake -DCMAKE_INSTALL_PREFIX=. .
     make install

   Detailed instructions can be found in this `Introduction to bv_maker`_ (login *brainvisa*, password *Soma2009*).

3. Create the configuration file for ``bv_maker`` at ``$HOME/.brainvisa/bv_maker.cfg``. Here is a minimal version of this file::

     [ source $HOME/brainvisa/source ]
       + brainvisa-cmake bug_fix
       + brainvisa-share bug_fix
       + soma-base bug_fix
       + soma-io bug_fix
       + aims-free bug_fix
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
- GSL_ (GNU Scientific Library).
- Boost_ version 1.49 or later.
- Python_ version 2.6 or later.
- CMake_ version 2.6 or later, with its extension ``brainvisa-cmake`` (distributed with BrainVISA_).
- Optional: the ``VipHomotopic`` command-line tool from the Morphologist image segmentation pipeline, distributed as a binary only tool with the BrainVISA_ installer.


Licence
=======

The source code of this work is placed under the CeCILL licence (see `<LICENCE.CeCILL.txt>`_). Compiled code that links to the GPL-licensed GSL_ forms a derivative work thereof, and thus must be redistributed under the GNU General Public Licence (see `<LICENCE.GPLv3.txt>`_).

.. _BrainVISA: http://brainvisa.info/
.. _GSL: http://www.gnu.org/software/gsl/
.. _Boost: http://www.boost.org/
.. _CMake: http://www.cmake.org/
.. _Python: https://www.python.org/
.. _Introduction to bv_maker: https://bioproj.extra.cea.fr/redmine/projects/brainvisa-devel/wiki/How_to_compile_BrainVISA_projects
.. _BrainVISA download page: http://brainvisa.info/web/download.html

.. Copyright Forschungszentrum Jülich GmbH (2016, 2017).
   Copyright Télécom ParisTech (2015, 2016).
   Copyright CEA (2014, 2015).
   Copyright Université Paris XI (2014).

   Author: Yann Leprince <yann.leprince@ylep.fr>.

   Copying and distribution of this file, with or without modification, are permitted in any medium without royalty provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.
