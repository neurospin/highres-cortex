# What is highres-cortex?

This is a collection of software designed to process high-resolution magnetic resonance images of the cerebral cortex.

If you use this work in an academic publication, **please cite** the relevant references, as listed in [`doc/references.bib`](doc/references.bib).


# Layout

The source is organized in three components, from high-level to low-level:

  - Scripts to do the high-level processing, under [`scripts/`](scripts/).
  - A Python package, under [`python/highres-cortex/`](python/highres-cortex/).
  - A library and compiled executables, written in C++, under [`src`](src/).


# Compiling

This package is compiled using the BrainVISA compilation system, which is an extension to CMake. Follow the [instructions at brainvisa.info](http://brainvisa.info/repository.html#use_brainvisa_sources) to clone the source tree. You need to include two lines in `bv_maker.cfg` in order to clone and compile `highres-cortex`. Here is a minimal `bv_maker.cfg`:

    [ source $HOME/brainvisa/source/trunk ]
      anatomist trunk
      git https://github.com/neurospin/highres-cortex.git master highres-cortex

    [ build $HOME/brainvisa/build/trunk ]
      anatomist trunk $HOME/brainvisa/source/trunk
      + $HOME/brainvisa/source/trunk/highres-cortex


# Dependencies

  - AIMS version 4.5 or later (in developement at the time of this writing, use the `trunk` version from SVN). This image processing library is distributed as part of [BrainVISA](http://brainvisa.info/).
  - [GSL (GNU Scientific Library)](http://www.gnu.org/software/gsl/).
  - [Boost](http://www.boost.org/) version 1.49 or later.
  - [CMake](http://www.cmake.org/) version 2.6 or later, with its extension `brainvisa-cmake` (distributed with [BrainVISA](http://brainvisa.info/)).
  - [Python](https://www.python.org/) version 2.6 or later.


# Licence

The source code of this work is placed under the [CeCILL licence](LICENCE.CeCILL.txt). Compiled code that links to the GPL-licenced [GSL (GNU Scientific Library)](http://www.gnu.org/software/gsl/) forms a derivated work thereof, and thus must be redistributed under the [GNU General Public Licence](LICENCE.GPLv3.txt).


# Legalese

Copyright CEA (2014).
Copyright Université Paris XI (2014).

Contributor: Yann Leprince <yann.leprince@ylep.fr>.

Copying and distribution of this file, with or without modification, are permitted in any medium without royalty provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.
