=========================================
Creating a Docker image with dependencies
=========================================

Such an image can be used to run tests in a controlled environment, or to create packages that are compatible with the official BrainVISA binary packages.

1. Download and install ``casa-distro`` according to its `documentation <http://brainvisa.info/casa-distro/>`.

2. Run ``casa_distro create_build_workflow distro=highres-cortex-test base_distro=opensource branch=latest_release system=ubuntu-12.04``

3. Optionally, customize configuration files in ``$HOME/casa_distro/highres-cortex-test/latest_release_ubuntu-12.04/conf``. Use the following ``bv_maker.cfg`` to include only minimal dependencies, which allows a faster build and a smaller resulting image::

    [ source $CASA_SRC ]
      brainvisa brainvisa-cmake $CASA_BRANCH
      brainvisa soma-base $CASA_BRANCH
      brainvisa soma-io $CASA_BRANCH
      brainvisa aims-free $CASA_BRANCH
      brainvisa morphologist-private $CASA_BRANCH
      brainvisa soma-workflow $CASA_BRANCH
      brainvisa capsul $CASA_BRANCH

    [ build $CASA_BUILD ]
      default_steps = configure build
      make_options = -j2
      build_type = RelWithDebInfo
      packaging_thirdparty = OFF
      clean_config = ON
      clean_build = ON
      brainvisa-cmake $CASA_BRANCH $CASA_SRC
      soma-base $CASA_BRANCH $CASA_SRC
      soma-io $CASA_BRANCH $CASA_SRC
      aims-free $CASA_BRANCH $CASA_SRC
      morphologist-private $CASA_BRANCH $CASA_SRC
      soma-workflow $CASA_BRANCH $CASA_SRC
      capsul $CASA_BRANCH $CASA_SRC

4. Run ``casa_distro bv_maker distro=highres-cortex-test branch=latest_release system=ubuntu-12.04``, which downloads and compiles the previously configured BrainVISA components in ``~/casa_distro/highres-cortex-test/latest_release_ubuntu-12.04``.

5. Install the software in a Docker container::

     docker run -v ~/casa_distro/highres-cortex-test/latest_release_ubuntu-12.04:/casa -it cati/casa-dev:ubuntu-12.04 /bin/bash

   Run the following commands in the Docker container:

   1. Install BrainVISA (under ``/usr/local``)::

        cd /casa/build
        make install-runtime install-dev

   2. BrainVISA installs Python packages in a non-standard path, so we need to link them so that Python can find them::

        cd /usr/local/lib/python2.7/dist-packages
        ln -s ../../../python/* .

6. Create a Docker image from the container::

     docker commit <container> casa-highres-cortex-test:ubuntu-12.04-latest_release
