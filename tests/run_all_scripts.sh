#! /bin/sh -e

: ${SCRIPTS_DIR:=$(dirname -- "$0")/../scripts}

run_script_in_subdir() {
    mkdir -p "$1"
    (
        cd "$1"
        "$SCRIPTS_DIR"/"$1"/"$2"
    )
}

run_script_in_subdir dist distmaps.sh
run_script_in_subdir heat heat.sh
run_script_in_subdir laplace-euclidean laplace-euclidean.sh
run_script_in_subdir upwind-euclidean upwind-euclidean.sh
run_script_in_subdir isovolume isovolume.sh
run_script_in_subdir column-regions column-regions.sh
