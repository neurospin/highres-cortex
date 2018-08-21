==========================
 Packaging highres-cortex
==========================

**TODO**: merge these instructions with *docker-test-image*.

Binary packages for highres-cortex are created with casa-distro. This procedure uses the singularity branch of casa-distro, which will become the default branch in the future.

This procedure is tested against BrainVISA release 4.6.1.

1. Initialize the build workflow (run the first time only)::

    casa_distro create distro_source=opensource distro_name=highres-cortex branch=latest_release system=ubuntu-12.04


2. Update ``$HOME/casa_distro/highres-cortex-test/latest_release_ubuntu-12.04/conf/bv_maker.cfg``::

    [ source $CASA_SRC ]
      brainvisa brainvisa-cmake $CASA_BRANCH
      brainvisa soma-base $CASA_BRANCH
    brainvisa soma-io $CASA_BRANCH
      brainvisa aims-free $CASA_BRANCH
      brainvisa soma-workflow $CASA_BRANCH
      brainvisa capsul $CASA_BRANCH
      brainvisa brainvisa-installer $CASA_BRANCH

    [ build $CASA_BUILD ]
    [ if os.getenv('WINEARCH') in ('win32', 'win64') ]
      cmake_options = -DBRAINVISA_CROSSCOMPILATION_DIR=$CASA_DEPS
    [ if os.getenv('WINEARCH') == 'win32' ]
      cross_compiling_prefix = i686-w64-mingw32
    [ else ]
    [if os.getenv('WINEARCH') == 'win64' ]
      cross_compiling_prefix = x86_64-w64-mingw32
    [ endif ]
    [ endif ]
    [ endif ]
      default_steps = configure build doc
      make_options = -j16
      build_type = Release
      packaging_thirdparty = ON
      clean_config = ON
      clean_build = ON
      test_ref_data_dir = $CASA_TESTS/ref
      test_run_data_dir = $CASA_TESTS/test
      brainvisa brainvisa-cmake $CASA_BRANCH $CASA_SRC
      brainvisa soma-base $CASA_BRANCH $CASA_SRC
      brainvisa soma-io $CASA_BRANCH $CASA_SRC
      brainvisa aims-free $CASA_BRANCH $CASA_SRC
      brainvisa soma-workflow $CASA_BRANCH $CASA_SRC
      brainvisa capsul $CASA_BRANCH $CASA_SRC
      brainvisa brainvisa-installer $CASA_BRANCH $CASA_SRC

    [ build $CASA_CUSTOM_BUILD ]
    [ if os.getenv('WINEARCH') in ('win32', 'win64') ]
      cmake_options = -DBRAINVISA_CROSSCOMPILATION_DIR=$CASA_DEPS
    [ if os.getenv('WINEARCH') == 'win32' ]
      cross_compiling_prefix = i686-w64-mingw32
    [ else ]
    [if os.getenv('WINEARCH') == 'win64' ]
      cross_compiling_prefix = x86_64-w64-mingw32
    [ endif ]
    [ endif ]
    [ endif ]
      make_options = -j16
      build_type = Release
      clean_config = ON
      clean_build = ON
      test_ref_data_dir = $CASA_TESTS/custom_ref
      test_run_data_dir = $CASA_TESTS/custom_test
      git https://github.com/neurospin/highres-cortex.git brainvisa_4.6.1 highres-cortex

3. Compile::

    casa_distro bv_maker distro=highres-cortex branch=latest_release system=ubuntu-12.04

4. Fix a conflicting ``sitecustomize`` directory (should be fixed in a future version of brainvisa-cmake)::

    casa_distro run distro=highres-cortex branch=latest_release system=ubuntu-12.04 rm -r /casa/custom/build/python/sitecustomize/

5. Run tests (``bv_maker test`` does not work for ``custom`` casa-distro projects in brainvisa-cmake 2.1.1)::

    casa_distro run distro=highres-cortex branch=latest_release system=ubuntu-12.04 /casa/build/bin/bv_env_test sh -e /casa/custom/build/build_files/highres-cortex/tests/test_env.sh python -m highres_cortex.test.test_capsul
    casa_distro run distro=highres-cortex branch=latest_release system=ubuntu-12.04 /casa/build/bin/bv_env_test sh -e /casa/custom/build/build_files/highres-cortex/tests/test_env.sh sh -e /casa/custom/src/highres-cortex/tests/test_all_scripts.sh

6. Create the package: **TODO**
