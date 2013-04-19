#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest
import numpy
import os

import lsst.utils.tests
import lsst.afw.geom
import lsst.afw.geom.ellipses as el
import lsst.shapelet.tractor
import lsst.shapelet.tests
import lsst.afw.image

# These parameters match those used to generate the check images; see
# tests/data/generate.py
GALAXY_RADIUS = 8.0
PSF_SIGMA = 2.0
E1 = 0.3
E2 = -0.2
PROFILES = [
    ("exp", 9, 8),
    ("ser2", 9, 8),
    ("ser3", 9, 8),
    ("dev", 9, 8),
    ]

CHECK_COMPONENT_IMAGES = False

class ProfileTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def testProfiles(self):
        psf0 = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE, PSF_SIGMA)
        psf0.getCoefficients()[:] = 1.0 / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
        psf = lsst.shapelet.MultiShapeletFunction()
        psf.getElements().push_back(psf0)
        psf.normalize()
        ellipse = el.Separable[el.Distortion, el.DeterminantRadius](E1, E2, GALAXY_RADIUS)
        for profile, nComponents, maxRadius in PROFILES:
            # check1 is the multi-Gaussian approximation, as convolved and evaluated by GalSim,
            check1 = lsst.afw.image.ImageD(os.path.join("tests", "data", profile + "-approx.fits")).getArray()
            xc = check1.shape[1] // 2
            yc = check1.shape[0] // 2
            xb = numpy.arange(check1.shape[1], dtype=float) - xc
            yb = numpy.arange(check1.shape[0], dtype=float) - yc
            xg, yg = numpy.meshgrid(xb, yb)

            basis = lsst.shapelet.tractor.loadBasis(profile, nComponents, maxRadius)
            builder = lsst.shapelet.MultiShapeletMatrixBuilderD(basis, psf, xg.ravel(), yg.ravel())
            image1 = numpy.zeros(check1.shape, dtype=float)
            matrix = image1.reshape(check1.size, 1)
            builder.build(matrix, el.Ellipse(ellipse))
            self.assertClose(check1, image1, plotOnFailure=True, rtol=5E-5, relTo=check1.max())
            msf = basis.makeFunction(el.Ellipse(ellipse, lsst.afw.geom.Point2D(xc, yc)),
                                     numpy.array([1.0], dtype=float))
            msf = msf.convolve(psf)
            image2 = numpy.zeros(check1.shape, dtype=float)
            msf.evaluate().addToImage(lsst.afw.image.ImageD(image2, False))
            self.assertClose(check1, image2, plotOnFailure=True, rtol=5E-5, relTo=check1.max())

            if profile == 'exp':
                # check2 is the exact profile, again by GalSim.
                # We only check exp against the exact profile.  The other approximations are less
                # accurate, and we only really need to test one.  The real measure of whether these
                # profiles are good enough is more complicated than what we can do in a unit test.
                check2 = lsst.afw.image.ImageD(
                    os.path.join("tests", "data", profile + "-exact.fits")
                ).getArray()
                self.assertClose(check2, image1, plotOnFailure=True, rtol=1E-3, relTo=check2.max())

            if CHECK_COMPONENT_IMAGES:
                # This was once useful for debugging test failures, and may be again, but it's
                # redundant with the above and requires putting more check images in git, so
                # it's disabled by default.
                for n, sf in enumerate(msf.getElements()):
                    check = lsst.afw.image.ImageD(
                        os.path.join("tests", "data", "%s-approx-%0d.fits" % (profile, n))
                    ).getArray()
                    image = numpy.zeros(check1.shape, dtype=float)
                    sf.evaluate().addToImage(lsst.afw.image.ImageD(image, False))
                    self.assertClose(check, image, plotOnFailure=True, rtol=5E-5, relTo=check1.max())


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ProfileTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
