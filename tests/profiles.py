#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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
    ("dev", 9, 8),
]

CHECK_COMPONENT_IMAGES = False


class ProfileTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def testRadii(self):
        """Check RadialProfile definitions of moments and half-light radii.
        """
        s = numpy.linspace(-20.0, 20.0, 1000)
        x, y = numpy.meshgrid(s, s)
        r = (x**2 + y**2)**0.5
        dxdy = (s[1] - s[0])**2
        for name in ["gaussian", "exp", "ser2", "luv", "lux"]:
            profile = lsst.shapelet.RadialProfile.get(name)
            z = profile.evaluate(r) * dxdy
            # lux and luv don't use the true half-light radius; instead they use the half-light radius
            # of the exp and dev profiles they approximate
            if not name.startswith("lu"):
                self.assertClose(z[r < 1].sum(), 0.5*z.sum(), rtol=0.01)
            # lhs of this comparison is the moments radius (using a sum approximation to the integral)
            self.assertClose(((z*x**2).sum() / z.sum())**0.5, profile.getMomentsRadiusFactor(), rtol=0.01)

    def testGaussian(self):
        """Test that the Gaussian profile's shapelet 'approximation' is actually exact.
        """
        profile = lsst.shapelet.RadialProfile.get("gaussian")
        r = numpy.linspace(0.0, 4.0, 100)
        z1 = profile.evaluate(r)
        basis = profile.getBasis(1)
        z2 = lsst.shapelet.tractor.evaluateRadial(basis, r, sbNormalize=True)[0, :]
        self.assertClose(z1, z2, rtol=1E-8)

    def testShapeletApproximations(self):
        psf0 = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE, PSF_SIGMA)
        psf0.getCoefficients()[:] = 1.0 / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
        psf = lsst.shapelet.MultiShapeletFunction()
        psf.getComponents().push_back(psf0)
        psf.normalize()
        ellipse = el.Separable[el.Distortion, el.DeterminantRadius](E1, E2, GALAXY_RADIUS)
        for name, nComponents, maxRadius in PROFILES:
            # check1 is the multi-Gaussian approximation, as convolved and evaluated by GalSim,
            check1 = lsst.afw.image.ImageD(os.path.join("tests", "data", name + "-approx.fits")).getArray()
            xc = check1.shape[1] // 2
            yc = check1.shape[0] // 2
            xb = numpy.arange(check1.shape[1], dtype=float) - xc
            yb = numpy.arange(check1.shape[0], dtype=float) - yc
            xg, yg = numpy.meshgrid(xb, yb)

            basis = lsst.shapelet.RadialProfile.get(name).getBasis(nComponents, maxRadius)
            builder = lsst.shapelet.MatrixBuilderD.Factory(xg.ravel(), yg.ravel(), basis, psf)()
            image1 = numpy.zeros(check1.shape, dtype=float)
            matrix = image1.reshape(check1.size, 1)
            builder(matrix, el.Ellipse(ellipse))
            self.assertClose(check1, image1, plotOnFailure=False, rtol=5E-5, relTo=check1.max())
            msf = basis.makeFunction(el.Ellipse(ellipse, lsst.afw.geom.Point2D(xc, yc)),
                                     numpy.array([1.0], dtype=float))
            msf = msf.convolve(psf)
            image2 = numpy.zeros(check1.shape, dtype=float)
            msf.evaluate().addToImage(lsst.afw.image.ImageD(image2, False))
            self.assertClose(check1, image2, plotOnFailure=False, rtol=5E-5, relTo=check1.max())

            if name == 'exp':
                # check2 is the exact profile, again by GalSim.
                # We only check exp against the exact profile.  The other approximations are less
                # accurate, and we only really need to test one.  The real measure of whether these
                # profiles are good enough is more complicated than what we can do in a unit test.
                check2 = lsst.afw.image.ImageD(
                    os.path.join("tests", "data", name + "-exact.fits")
                ).getArray()
                self.assertClose(check2, image1, plotOnFailure=False, rtol=1E-3, relTo=check2.max())

            if CHECK_COMPONENT_IMAGES:
                # This was once useful for debugging test failures, and may be again, but it's
                # redundant with the above and requires putting more check images in git, so
                # it's disabled by default.
                for n, sf in enumerate(msf.getComponents()):
                    check = lsst.afw.image.ImageD(
                        os.path.join("tests", "data", "%s-approx-%0d.fits" % (name, n))
                    ).getArray()
                    image = numpy.zeros(check1.shape, dtype=float)
                    sf.evaluate().addToImage(lsst.afw.image.ImageD(image, False))
                    self.assertClose(check, image, plotOnFailure=False, rtol=5E-5, relTo=check1.max())


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
