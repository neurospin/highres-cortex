#! /bin/sh

# This script will download highres-cortex and its BrainVISA dependencies, and
# try to build them from source. It will create a fresh directory under /tmp
# and do everything under this directory. The directory will *not* be deleted
# by the script, so that you can try to run the compiled programs by yourself.

die() {
    echo "Fatal error: $1" >&2
    echo "Test directory $TEST_DIR left behind." >&2
    exit 1
}

TEST_DIR=$(mktemp -d) || die "cannot create temporary directory"

echo "Testing checkout and build under $TEST_DIR"

BRAINVISA_CMAKE_SVN=https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/branches/bug_fix

svn export --username brainvisa --password Soma2009 "$BRAINVISA_CMAKE_SVN" $TEST_DIR/brainvisa-cmake || die "cannot fetch brainvisa-cmake"

(
    set -e
    cd $TEST_DIR/brainvisa-cmake
    cmake -DCMAKE_INSTALL_PREFIX=. .
    make install
) || die "cannot bootstrap brainvisa-cmake"

cat <<EOF > $TEST_DIR/bv_maker.cfg
[ source $TEST_DIR/source ]
  + brainvisa-cmake bug_fix
  + soma-base bug_fix
  + soma-io bug_fix
  + aims-free bug_fix
  git https://github.com/neurospin/highres-cortex.git master highres-cortex

[ build $TEST_DIR/build ]
  build_type = Release
  make_options = -j$(nproc 2>/dev/null || echo 1)
  brainvisa-cmake bug_fix $TEST_DIR/source
  soma-base bug_fix $TEST_DIR/source
  soma-io bug_fix $TEST_DIR/source
  aims-free bug_fix $TEST_DIR/source
  + $TEST_DIR/source/highres-cortex
EOF


BV_MAKER="$TEST_DIR/brainvisa-cmake/bin/bv_maker -c $TEST_DIR/bv_maker.cfg"

$BV_MAKER sources || die "fetching the sources failed"
$BV_MAKER configure || die "configuring with CMake failed"
$BV_MAKER build || die "build failed"
$BV_MAKER test || die "test failed"

echo "Success."
echo "Test directory $TEST_DIR left behind."
