# -*- coding: utf-8 -*-
# Copyright CEA (2021).
# Copyright Forschungszentrum Jülich GmbH (2017).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# This file is part of highres-cortex, a collection of software designed
# to process high-resolution magnetic resonance images of the cerebral
# cortex.
#
# This software is governed by the CeCILL licence under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# licence as circulated by CEA, CNRS and INRIA at the following URL:
# <http://www.cecill.info/>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the licence, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of scientific
# software, that may mean that it is complicated to manipulate, and that
# also therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL licence and that you accept its terms.

from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile
import unittest

import capsul.api

from . import synthetic_data
from . import compare_with_reference


class SphereTestCase(unittest.TestCase):
    def setUp(self):
        try:
            self.test_dir = tempfile.mkdtemp(
                prefix="highres-cortex-capsul-tests")
            synthetic_data.write_sphere_and_reference_result(
                1, 4, 0.3, dir=self.test_dir)

            self.result_comp = compare_with_reference.ResultComparator(
                self.test_dir)

            self.capsul_engine = capsul.api.Capsul().engine()

            p1 = capsul.api.executable(
                "highres_cortex.capsul.processes.BinarizeCortex")
            p1.classif = os.path.join(self.test_dir, "classif.nii.gz")
            p1.output_image = os.path.join(self.test_dir, "cortex_mask.nii.gz")
            with self.capsul_engine as ce:
                ce.run(p1)
        except BaseException:
            if hasattr(self, "test_dir"):
                shutil.rmtree(self.test_dir)
            raise
        if os.environ.get('KEEP_TEMPORARY'):
            print('highres-cortex test directory is {0}'.format(self.test_dir))

    def tearDown(self):
        if not os.environ.get('KEEP_TEMPORARY'):
            shutil.rmtree(self.test_dir)

    def test_laplacian(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.processes.Laplacian")
        p.classif = os.path.join(self.test_dir, "classif.nii.gz")
        p.precision = 0.001
        p.typical_cortical_thickness = 3
        p.laplace_field = os.path.join(self.test_dir, "laplacian.nii.gz")
        with self.capsul_engine as ce:
            ce.run(p)
        res = self.result_comp.ensure_max_rms_error(
            "laplacian.nii.gz", 0.017,
            reference_file="reference_laplacian.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_filtered_sumcurvs(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.filtered_sumcurvs")
        p.input = os.path.join(self.test_dir, "reference_laplacian.nii.gz")
        p.mode = "sum"
        p.output = os.path.join(self.test_dir, "curvature.nii.gz")
        with self.capsul_engine as ce:
            ce.run(p)
        res = self.result_comp.ensure_max_rms_error(
            "curvature.nii.gz", 0.067,
            reference_file="reference_curvature.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_advect_euclidean(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.processes.EuclideanAdvectionAlongGradient")
        p.domain = os.path.join(self.test_dir, "cortex_mask.nii.gz")
        p.grad_field = os.path.join(
            self.test_dir, "reference_laplacian.nii.gz")
        p.step_size = 0.05
        p.upfield = False
        p.output_length = os.path.join(
            self.test_dir, "euclidean_adv_toward_white.nii.gz")
        p()
        res = self.result_comp.ensure_max_rms_error(
            "euclidean_adv_toward_white.nii.gz", 0.075,
            reference_file="reference_distwhite.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_upwind_euclidean(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.processes.EuclideanUpwindingAlongGradient")
        p.domain = os.path.join(self.test_dir, "classif.nii.gz")
        p.field = os.path.join(
            self.test_dir, "reference_laplacian.nii.gz")
        p.downfield = True
        p.origin_label = 200
        p.output = os.path.join(
            self.test_dir, "euclidean_upw_toward_white.nii.gz")
        p()
        res = self.result_comp.ensure_max_rms_error(
            "euclidean_upw_toward_white.nii.gz", 0.22,
            reference_file="reference_distwhite.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_equivolumetric_pipeline(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.isovolume")
        p.classif = os.path.join(self.test_dir, "classif.nii.gz")
        p.advection_step_size = 0.05
        p.equivolumetric_depth = os.path.join(
            self.test_dir, "equivolumetric_depth.nii.gz")
        p()
        res = self.result_comp.ensure_max_rms_error(
            "equivolumetric_depth.nii.gz", 0.028,
            reference_file="reference_equivolumic.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_thickness_adv_pipeline(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.thickness_adv")
        p.classif = os.path.join(self.test_dir, "classif.nii.gz")
        p.advection_step_size = 0.05
        p.thickness_image = os.path.join(
            self.test_dir, "thickness_adv.nii.gz")
        p.equidistant_depth = os.path.join(
            self.test_dir, "equidistant_depth_adv.nii.gz")
        p()
        res = self.result_comp.ensure_max_rms_error(
            "thickness_adv.nii.gz", 0.12,
            reference_file="reference_thickness.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")
        res = self.result_comp.ensure_max_rms_error(
            "equidistant_depth_adv.nii.gz", 0.017,
            reference_file="reference_euclidean.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_thickness_upw_pipeline(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.thickness_upw")
        p.classif = os.path.join(self.test_dir, "classif.nii.gz")
        p.thickness_image = os.path.join(
            self.test_dir, "thickness_upw.nii.gz")
        p.equidistant_depth = os.path.join(
            self.test_dir, "equidistant_depth_upw.nii.gz")
        p()
        res = self.result_comp.ensure_max_rms_error(
            "thickness_upw.nii.gz", 0.27,
            reference_file="reference_thickness.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")
        res = self.result_comp.ensure_max_rms_error(
            "equidistant_depth_upw.nii.gz", 0.024,
            reference_file="reference_euclidean.nii.gz")
        self.assertTrue(res, msg="RMS error is too high")

    def test_traverses_pipeline(self):
        p = capsul.api.executable(
            "highres_cortex.capsul.traverses")
        p.classif = os.path.join(self.test_dir, "classif.nii.gz")
        p.cortical_traverses = os.path.join(
            self.test_dir, "traverses.nii.gz")
        p()


if __name__ == "__main__":
    unittest.main()
