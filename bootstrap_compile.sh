#! /bin/sh


echo "This script will download and build highres-cortex and its dependencies."


die() {
    echo "Fatal error: $1" >&2
    exit 1
}


#####################################
# Check if dependencies are installed
#####################################
if lsb_release -c | grep trusty 2>&1 >/dev/null; then
    # (Ubuntu 14.04 "trusty")
    packages="subversion git cmake make gcc g++ pkg-config libsigc++-2.0-dev libxml2-dev python-sip-dev libboost-dev libgsl0-dev python-numpy python2.7-dev"
    run_aptget=false

    for package in $PACKAGES; do
        if dpkg -s $package 2>&1 >/dev/null; then
            run_aptget=true
            break
        fi
    done

    if $run_aptget; then
        sudo apt-get --no-install-recommends install $packages || die "cannot carry on due to missing dependencies"
    fi
else
    # Unknown distribution/version
    echo "Your distribution is unknown. Please check README.rst"
    echo "and install the required dependencies."
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

echo "Will now bootstrap highres-cortex under $base_dir"


###########################
# Bootstrap brainvisa-cmake
###########################

BRAINVISA_CMAKE_SVN=https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/branches/bug_fix

svn export --non-interactive --username brainvisa --password Soma2009 "$BRAINVISA_CMAKE_SVN" $base_dir/brainvisa-cmake || die "cannot download brainvisa-cmake"

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
EOF


BV_MAKER="$base_dir/brainvisa-cmake/bin/bv_maker"

$BV_MAKER sources || die "fetching the sources failed"
$BV_MAKER configure || die "configuring with CMake failed"
$BV_MAKER build || die "build failed"

echo
echo "You can now use highres-cortex."
echo "Remember to set up the environment first with:"
echo "  . $base_dir/build/bin/bv_env.sh $base_dir/build"
echo
echo "You can update highres-cortex by running:"
echo "  $base_dir/build/bin/bv_maker"
echo
