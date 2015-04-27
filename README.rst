================
 highres-cortex
================

This is a collection of software designed to process high-resolution magnetic resonance images of the cerebral cortex.

If you use this work in an academic publication, **please cite** the relevant references (see `<doc/references.bib>`_).


Layout
------

The source is organized in three components, from high-level to low-level:

- Scripts to do the high-level processing, under `<scripts/>`_.
- A Python package, under `<python/highres_cortex/>`_.
- A library and compiled executables, written in C++, under `<src/>`_.


Compiling
---------

This package is compiled using the BrainVISA_ compilation system, which is an extension to CMake_. Follow the instructions at http://brainvisa.info/repository.html#use_brainvisa_sources to clone the source tree. You need to include two lines in ``bv_maker.cfg`` in order to clone and compile ``highres-cortex``. Here is a minimal ``bv_maker.cfg``::

    [ source $HOME/brainvisa/source/trunk ]
      + anatomist trunk
      git https://github.com/neurospin/highres-cortex.git master highres-cortex

    [ build $HOME/brainvisa/build/trunk ]
      anatomist trunk $HOME/brainvisa/source/trunk
      + $HOME/brainvisa/source/trunk/highres-cortex


Dependencies
------------

- AIMS version 4.5 or later (in development at the time of this writing, use the ``trunk`` version from SVN). This image processing library is distributed as part of BrainVISA_.
- The ``VipHomotopic`` tool from the Morphologist image segmentation pipeline, distributed as part of BrainVISA_.
- GSL_ (GNU Scientific Library).
- Boost_ version 1.49 or later.
- CMake_ version 2.6 or later, with its extension ``brainvisa-cmake`` (distributed with BrainVISA_).
- Python_ version 2.6 or later.


Licence
-------

The source code of this work is placed under the CeCILL licence (see `<LICENCE.CeCILL.txt>`_). Compiled code that links to the GPL-licensed GSL_ forms a derivative work thereof, and thus must be redistributed under the GNU General Public Licence (see `<LICENCE.GPLv3.txt>`_).


.. Copyright CEA (2014).
   Copyright Université Paris XI (2014).

   Contributor: Yann Leprince <yann.leprince@ylep.fr>.

   Copying and distribution of this file, with or without modification, are permitted in any medium without royalty provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.

.. _BrainVISA: http://brainvisa.info/
.. _GSL: http://www.gnu.org/software/gsl/
.. _Boost: http://www.boost.org/
.. _CMake: http://www.cmake.org/
.. _Python: https://www.python.org/
