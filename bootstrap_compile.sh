#! /bin/sh


# Helper functions
die() {
    echo "Fatal error: $1" >&2
    exit 1
}

msg() {
    tput smso  # bold
    echo "$@"
    tput rmso  # bold off
}


msg "This script will download and build highres-cortex and its dependencies."


#####################################
# Check if dependencies are installed
#####################################
if lsb_release -c | grep trusty 2>&1 >/dev/null; then
    # (Ubuntu 14.04 "trusty")
    # libqt4-sql-sqlite is not needed, but warnings show up when it is not present
    packages="subversion ca-certificates git cmake make gcc g++ gfortran pkg-config libblitz0-dev libsigc++-2.0-dev libxml2-dev libqt4-dev libboost-dev zlib1g-dev libtiff-dev libgsl0-dev python2.7-dev python-sip-dev python-numpy python-six libqt4-sql-sqlite"
    run_aptget=false

    for package in $packages; do
        if ! dpkg -s $package 2>&1 >/dev/null; then
            run_aptget=true
            break
        fi
    done

    if [ "$run_aptget" = true ]; then
        sudo apt-get --no-install-recommends install $packages || die "cannot carry on due to missing dependencies"
    fi
else
    # Unknown distribution/version
    msg "Your distribution is unknown. Please check README.rst"
    msg "and install the required dependencies yourself."
fi



########################################
# Ask for base directory for compilation
########################################

# default value
base_dir=$HOME/brainvisa

printf 'Where do you want to bootstrap highres-cortex? [Default: %s] ' "$HOME/brainvisa"
read base_dir

if [ -z "$base_dir" ]; then
    base_dir=$HOME/brainvisa
fi

mkdir -p -- "$base_dir" || die "cannot create directory '$base_dir'"
cd -- "$base_dir" || die "cannot enter directory '$base_dir'"

msg "Will now bootstrap highres-cortex under $base_dir"


###########################
# Bootstrap brainvisa-cmake
###########################

brainvisa_cmake_svn=https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/branches/bug_fix

svn export --non-interactive --username brainvisa --password Soma2009 "$brainvisa_cmake_svn" $base_dir/brainvisa-cmake || die "cannot download brainvisa-cmake"

(
    set -e
    cd $base_dir/brainvisa-cmake
    cmake -DCMAKE_INSTALL_PREFIX=. .
    make install
) || die "cannot bootstrap brainvisa-cmake"


############################################################
# Use brainvisa-cmake to download and compile highres-cortex
############################################################

mkdir -p "$HOME/.brainvisa"
cat <<EOF > "$HOME/.brainvisa/bv_maker.cfg"
[ source $base_dir/source ]
  + brainvisa-cmake bug_fix
  + brainvisa-share bug_fix
  + soma-base bug_fix
  + soma-io bug_fix
  + aims-free bug_fix
  git https://github.com/neurospin/highres-cortex.git master highres-cortex

[ build $base_dir/build ]
  build_type = Release
  make_options = -j$(nproc 2>/dev/null || echo 1)
  brainvisa-cmake bug_fix $base_dir/source
  brainvisa-share bug_fix $base_dir/source
  soma-base bug_fix $base_dir/source
  soma-io bug_fix $base_dir/source
  aims-free bug_fix $base_dir/source
  + $base_dir/source/highres-cortex
EOF


BV_MAKER=$base_dir/brainvisa-cmake/bin/bv_maker

msg
msg "******************************************************"
msg "* The next step (fetching sources of highres-cortex  *"
msg "* dependencies) may ask you for identification. Use: *"
msg "*    Username: brainvisa                             *"
msg "*    Password: Soma2009                              *"
msg "******************************************************"
msg

"$BV_MAKER" --username brainvisa sources || die "fetching the sources failed"
"$BV_MAKER" configure || die "configuring with CMake failed"
"$BV_MAKER" build || die "build failed"

# This is not needed anymore, use $base_dir/build/bin/bv_maker instead
rm -rf "$base_dir"/brainvisa-cmake

msg
msg "You can now use highres-cortex."
msg "Remember to set up the environment first with:"
msg "  . $base_dir/build/bin/bv_env.sh $base_dir/build"
msg
msg "You can update highres-cortex by running:"
msg "  $base_dir/build/bin/bv_maker"
msg
