#! /bin/sh

# This script will run interactively if run from a terminal. Otherwise, it
# should run non-interactively and use the default settings (e.g. install to
# $HOME/brainvisa, do not customize bv_maker.cfg).

# Helper functions
die() {
    msg "Fatal error: $1"
    exit 1
}
bold=$(tput smso) offbold=$(tput rmso) 2>/dev/null
msg() {
    echo "${bold}$*${offbold}" >&2
}
sudo_root() {
    # Run the command with sudo if available
    if type sudo >/dev/null 2>&1; then
        sudo -- "$@"
    else
        "$@"
    fi
}


msg "This script will download and build highres-cortex and its dependencies."


#####################################
# Check if dependencies are installed
#####################################

# Get variables identifying distribution
ID= ID_LIKE= PRETTY_NAME= VERSION_ID=
. /etc/os-release >/dev/null 2>&1 || . /usr/lib/os-release >/dev/null 2>&1

packages_installed=false
if [ "$ID" = ubuntu ] || [ "$ID_LIKE" = ubuntu ]; then
    packages=
    if [ "$VERSION_ID" = 16.04 ]; then
        packages="subversion ca-certificates git cmake make gcc g++ gfortran \
                  gcc-4.9 g++-4.9 gfortran-4.9 pkg-config libblitz0-dev \
                  libsigc++-2.0-dev libxml2-dev libqt4-dev libboost-dev \
                  zlib1g-dev libgsl-dev python2.7-dev python-sip-dev \
                  python-numpy python-six libqt4-sql-sqlite"
    elif [ "$VERSION_ID" = 14.04 ]; then
        # libqt4-sql-sqlite is not needed, but warnings show up when it is not
        # present
        packages="subversion ca-certificates git cmake make gcc g++ gfortran \
                  pkg-config libblitz0-dev libsigc++-2.0-dev libxml2-dev \
                  libqt4-dev libboost-dev zlib1g-dev libgsl0-dev \
                  python2.7-dev python-sip-dev python-numpy python-six \
                  libqt4-sql-sqlite"
    fi

    if [ -n "$packages" ]; then
        if ! dpkg-query -s $packages >/dev/null 2>&1; then
            msg "Will now install missing dependencies using apt-get"
            opt=
            if ! [ -t 0 ]; then
                # if running non-interactively (not attached to a terminal)
                opt="-y -q"
            fi
            sudo_root apt-get $opt --no-install-recommends install $packages \
                || die "missing dependencies could not be installed"
            packages_installed=true
        fi
    fi
fi
if ! [ $packages_installed = true ]; then
    # Unknown distribution/version
    msg "Your distribution ($PRETTY_NAME) is unknown."
    msg "You will have to install the dependencies yourself (see README.rst)"
fi


########################################
# Ask for base directory for compilation
########################################

# default value
base_dir=$HOME/brainvisa

printf 'Where do you want to bootstrap highres-cortex? [Default: %s] ' \
       "$HOME/brainvisa"
read base_dir

if [ -z "$base_dir" ]; then
    base_dir=$HOME/brainvisa
fi

mkdir -p -- "$base_dir" || die "cannot create directory '$base_dir'"
cd -- "$base_dir" || die "cannot enter directory '$base_dir'"

msg "Will now bootstrap highres-cortex under $PWD"


###########################
# Bootstrap brainvisa-cmake
###########################

if [ -x "$base_dir/brainvisa-cmake/bin/bv_maker" ]; then
    bv_maker=$base_dir/brainvisa-cmake/bin/bv_maker
    msg "
An existing checkout of brainvisa-cmake was detected at
$base_dir/brainvisa-cmake/bin/bv_maker

This script will use it to continue. If this causes problems, you can
delete the $base_dir/brainvisa-cmake directory
and run this script again to make a clean checkout."
else
    brainvisa_cmake_svn=https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/branches/bug_fix

    svn export --non-interactive --username brainvisa --password Soma2009 "$brainvisa_cmake_svn" $base_dir/brainvisa-cmake || die "cannot download brainvisa-cmake"

    (set -e
        cd "$base_dir"/brainvisa-cmake
        cmake -DCMAKE_INSTALL_PREFIX=. .
        make install
    ) || die "cannot bootstrap brainvisa-cmake"

    bv_maker=$base_dir/brainvisa-cmake/bin/bv_maker
fi


############################################################
# Use brainvisa-cmake to download and compile highres-cortex
############################################################

bv_maker_cfg_contents="\
[ source $base_dir/source ]
  + brainvisa-cmake bug_fix
  + soma-base bug_fix
  + soma-io bug_fix
  + aims-free bug_fix
  git https://github.com/neurospin/highres-cortex.git master highres-cortex

[ build $base_dir/build ]
  build_type = Release
  make_options = -j$(nproc 2>/dev/null || echo 1)
  brainvisa-cmake bug_fix $base_dir/source
  soma-base bug_fix $base_dir/source
  soma-io bug_fix $base_dir/source
  aims-free bug_fix $base_dir/source
  + $base_dir/source/highres-cortex
"

mkdir -p "$HOME/.brainvisa"
set -C  # do not overwrite existing configuration

echo "$bv_maker_cfg_contents" > "$HOME/.brainvisa/bv_maker.cfg" 2>/dev/null
if [ $? -ne 0 ] ; then
    msg "\
A bv_maker configuration already exists at $HOME/.brainvisa/bv_maker.cfg

Make sure that it contains the following:

$bv_maker_cfg_contents
"
fi

msg "\
Please adapt $HOME/.brainvisa/bv_maker.cfg now if needed, then
press Enter to continue"
read dummy


msg "
******************************************************
* The next step (fetching sources of highres-cortex  *
* dependencies) may ask you for identification. Use: *
*    Username: brainvisa                             *
*    Password: Soma2009                              *
******************************************************
"

"$bv_maker" --username brainvisa sources || die "fetching the sources failed"
"$bv_maker" configure || die "configuring with CMake failed"
"$bv_maker" build || die "build failed"

# This is not needed anymore, use $base_dir/build/bin/bv_maker instead
rm -rf "$base_dir"/brainvisa-cmake

msg "
You can now use highres-cortex.
Remember to set up the environment first with:
  . $base_dir/build/bin/bv_env.sh $base_dir/build

You can update highres-cortex by running:
  $base_dir/build/bin/bv_maker
"

exit 0
