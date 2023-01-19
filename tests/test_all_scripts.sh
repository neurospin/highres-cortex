#! /bin/sh -e
#
# Copyright Forschungszentrum Jülich GmbH (2017).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved. This file is offered as-is,
# without any warranty.

test_dir=$(dirname -- "$0")
[ -f "$test_dir/run_all_scripts.sh" ] || {
    echo "$0: cannot find the directory containing run_all_scripts.sh" >&2
    exit 1
}

# Create a temporary directory, make sure that it will be deleted on exit
tmpdir=
cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT
trap 'cleanup; trap - HUP EXIT; kill -HUP $$' HUP
trap 'cleanup; trap - INT EXIT; kill -INT $$' INT
trap 'cleanup; trap - TERM EXIT; kill -TERM $$' TERM
trap 'trap - QUIT EXIT; kill -QUIT $$' QUIT
tmpdir=$(mktemp -d -t hrcortex_XXXX)

# The Python interpreter below is called through bv_env in order to
# restore the values of the environment variables DYLD_LIBRARY_PATH
# and friends on Mac OS 10.11+, which are not inherited by the shell
# interpreter because it is protected (it resides in the protected
# directory /bin).

cd -- "$tmpdir"
bv_env python -m highres_cortex.test.synthetic_data 5 3 0.3
"$test_dir"/run_all_scripts.sh

bv_env python - <<EOF
import sys
from highres_cortex.test.compare_with_reference import *
c = ResultComparator(".")
success = c.ensure_max_rms_errors([
    ("heat/heat.nii.gz", 0.021),
    ("isovolume/equivolumic_depth.nii.gz", 0.020),
    ("laplace-euclidean/pial-fraction.nii.gz", 0.017),
    ("laplace-euclidean/total-length.nii.gz", 0.150),
    ("upwind-euclidean/pial-fraction.nii.gz", 0.023),
    ("upwind-euclidean/total-length.nii.gz", 0.190),
])
sys.exit(0 if success else 1)
EOF
