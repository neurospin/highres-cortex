#! /bin/sh -e

python make_test_classif.py


SCRIPTS_DIR=../scripts

cp -R "$SCRIPTS_DIR"/* .

(cd dist && ./distmaps.sh)
(cd heat && ./heat.sh)
(cd laplace-euclidean && ./laplace-euclidean.sh)
(cd upwind-euclidean && ./upwind-euclidean.sh)
(cd isovolume && ./isovolume.sh)
(cd column-regions && ./column-regions.sh)
