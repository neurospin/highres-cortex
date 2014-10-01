#! /bin/sh

die() {
    echo "Fatal error: $1" >&2
    echo "Test directory $TEST_DIR left behind." >&2
    exit 1
}

TEST_DIR=$(mktemp -d) || die "cannot create temporary directory"

echo "Testing checkout and build under $TEST_DIR"

BRAINVISA_CMAKE_SVN=https://bioproj.extra.cea.fr/neurosvn/brainvisa/development/brainvisa-cmake/trunk

svn export --username brainvisa --password Soma2009 "$BRAINVISA_CMAKE_SVN" $TEST_DIR/brainvisa-cmake || die "cannot fetch brainvisa-cmake"

(
    set -e
    cd $TEST_DIR/brainvisa-cmake
    cmake -DCMAKE_INSTALL_PREFIX=. .
    make install
) || die "cannot build brainvisa-cmake"

cat <<EOF > $TEST_DIR/bv_maker.cfg
[ source $TEST_DIR/source ]
  + anatomist trunk
  - brainvisa-share
  git https://github.com/neurospin/highres-cortex.git master highres-cortex

[ build $TEST_DIR/build ]
  anatomist trunk $TEST_DIR/source
  + $TEST_DIR/source/highres-cortex
EOF


BV_MAKER="$TEST_DIR/brainvisa-cmake/bin/bv_maker -c $TEST_DIR/bv_maker.cfg"

$BV_MAKER sources || die "fetching the sources failed"
$BV_MAKER configure || die "configuring with CMake failed"
$BV_MAKER build || die "build failed"

echo "Success."
echo "Test directory $TEST_DIR left behind."
