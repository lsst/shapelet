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

from __future__ import print_function
import unittest
import numpy
import cPickle

import lsst.utils.tests
import lsst.afw.geom.ellipses as ellipses
import lsst.shapelet.tests
import lsst.afw.image

numpy.random.seed(500)


class MultiShapeletTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def testMoments(self):
        x = numpy.linspace(-50, 50, 1001)
        y = x
        function = self.makeRandomMultiShapeletFunction()
        x = numpy.linspace(-10, 10, 101)
        y = x
        z = self.makeImage(function, x, y)
        self.checkMoments(function, x, y, z)

    def testPickle(self):
        function1 = self.makeRandomMultiShapeletFunction()
        s = cPickle.dumps(function1, protocol=2)
        function2 = cPickle.loads(s)
        for component1, component2 in zip(function1.getComponents(), function2.getComponents()):
            self.assertEqual(component1.getOrder(), component2.getOrder())
            self.assertEqual(component1.getBasisType(), component2.getBasisType())
            self.assertClose(component1.getEllipse().getParameterVector(),
                             component2.getEllipse().getParameterVector())
            self.assertClose(component1.getCoefficients(), component2.getCoefficients())

    def testConvolveGaussians(self):
        sigma1 = [lsst.afw.geom.ellipses.Quadrupole(6.0, 5.0, 2.0),
                  lsst.afw.geom.ellipses.Quadrupole(8.0, 10.0, -1.0)]
        sigma2 = [lsst.afw.geom.ellipses.Quadrupole(7.0, 12.0, -2.0),
                  lsst.afw.geom.ellipses.Quadrupole(7.0, 9.0, 1.0)]
        alpha1 = [0.6, 0.4]
        alpha2 = [0.35, 0.65]
        sigma3 = []
        alpha3 = []

        def makeMultiShapeletFunction(alpha, sigma):
            msf = lsst.shapelet.MultiShapeletFunction()
            for a, s in zip(alpha, sigma):
                f = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE,
                                                   lsst.afw.geom.ellipses.Ellipse(s))
                f.getCoefficients()[0] = a / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
                msf.getComponents().push_back(f)
            return msf
        for a1, s1 in zip(alpha1, sigma1):
            for a2, s2 in zip(alpha2, sigma2):
                sigma3.append(lsst.afw.geom.ellipses.Quadrupole(s1.getIxx() + s2.getIxx(),
                                                                s1.getIyy() + s2.getIyy(),
                                                                s1.getIxy() + s2.getIxy()))
                alpha3.append(a1 * a2)
        msf1 = makeMultiShapeletFunction(alpha1, sigma1)
        msf2 = makeMultiShapeletFunction(alpha2, sigma2)
        msf3a = makeMultiShapeletFunction(alpha3, sigma3)
        msf3b = msf1.convolve(msf2)

        # Compare the parameters of the MultiShapeletFunctions
        self.compareMultiShapeletFunctions(msf3a, msf3b)

        # Just to be extra sure, we test that the images are also the same
        x = numpy.linspace(-20, 20, 41)
        y = numpy.linspace(-20, 20, 41)
        image3a = self.makeImage(msf3a, x, y)
        image3b = self.makeImage(msf3b, x, y)
        self.assertClose(image3a, image3b)

        # Now we test against two test images: one implemented right here with numpy calls...
        xg, yg = numpy.meshgrid(x, y)

        def evalMultiGaussian(alpha, sigma):
            matQ = numpy.array([[sigma.getIxx(), sigma.getIxy()],
                                [sigma.getIxy(), sigma.getIyy()]],
                               dtype=float)
            invQ = numpy.linalg.inv(matQ)
            norm = alpha / numpy.linalg.det(2.0 * numpy.pi * matQ)**0.5
            return norm * numpy.exp(-0.5 * (invQ[0, 0]*xg**2 + 2.0*invQ[0, 1]*xg*yg + invQ[1, 1]*yg**2))
        image3c = numpy.zeros(xg.shape, dtype=float)
        for a, s in zip(alpha3, sigma3):
            image3c += evalMultiGaussian(a, s)
        self.assertClose(image3c, image3a, rtol=1E-6, relTo=numpy.max(image3c),
                         printFailures=True, plotOnFailure=False)

        # And the second produced by GalSim
        if False:
            # Print inputs to screen so we can make a test image with GalSim (see tests/data/generate.py)
            # Output can be pasted into that file to generate the check image.
            def printForGalSim(alpha, sigma):
                print("galsim.Add([")
                for a, s in zip(alpha, sigma):
                    e = lsst.afw.geom.ellipses.Separable[lsst.afw.geom.ellipses.Distortion,
                                                         lsst.afw.geom.ellipses.DeterminantRadius](s)
                    print("    makeGaussian(flux=%f, e1=%8.8f, e2=%8.8f, sigma=%8.8f),"
                          % (a, e.getE1(), e.getE2(), e.getRadius()))
                print("])")
            printForGalSim(alpha1, sigma1)
            printForGalSim(alpha2, sigma2)
        image3d = lsst.afw.image.ImageF("tests/data/gaussians.fits").getArray().astype(float)
        self.assertClose(image3d, image3a, rtol=1E-6, relTo=numpy.max(image3d),
                         printFailures=True, plotOnFailure=False)

    def testBasisNormalize(self):
        def makePositiveMatrix(*shape):
            """Return a random basis matrix, but with a lot of power
            in the zeroth component to ensure the integral is positve."""
            a = numpy.random.randn(*shape)
            a[0, :] += 4.0
            return a
        basis = lsst.shapelet.MultiShapeletBasis(2)
        basis.addComponent(0.5, 1, makePositiveMatrix(3, 2))
        basis.addComponent(1.0, 2, makePositiveMatrix(6, 2))
        basis.addComponent(1.2, 0, makePositiveMatrix(1, 2))
        basis.normalize()
        for n in range(2):
            coefficients = numpy.zeros(2, dtype=float)
            coefficients[n] = 1.0
            msf = basis.makeFunction(lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes()),
                                     coefficients)
            self.assertClose(msf.evaluate().integrate(), 1.0)

    def testBasisScale(self):
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes())
        basis = lsst.shapelet.MultiShapeletBasis(2)
        basis.addComponent(0.5, 1, numpy.random.randn(3, 2))
        basis.addComponent(1.0, 2, numpy.random.randn(6, 2))
        basis.addComponent(1.2, 0, numpy.random.randn(1, 2))
        msf1 = [basis.makeFunction(ellipse, self.makeUnitVector(i, 2)) for i in range(2)]
        basis.scale(2.0)
        ellipse.getCore().scale(0.5)
        msf2 = [basis.makeFunction(ellipse, self.makeUnitVector(i, 2)) for i in range(2)]
        for a, b in zip(msf1, msf2):
            self.compareMultiShapeletFunctions(a, b)

    def testBasisMerge(self):
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes())
        basis1 = lsst.shapelet.MultiShapeletBasis(2)
        basis1.addComponent(0.5, 1, numpy.random.randn(3, 2))
        basis1.addComponent(1.0, 2, numpy.random.randn(6, 2))
        basis1.addComponent(1.2, 0, numpy.random.randn(1, 2))
        basis2 = lsst.shapelet.MultiShapeletBasis(3)
        basis2.addComponent(0.4, 1, numpy.random.randn(3, 3))
        basis2.addComponent(1.1, 2, numpy.random.randn(6, 3))
        basis2.addComponent(1.6, 0, numpy.random.randn(1, 3))
        basis3 = lsst.shapelet.MultiShapeletBasis(basis1)
        basis3.merge(basis2)
        self.assertEqual(basis3.getSize(), 5)
        msf1 = [basis1.makeFunction(ellipse, self.makeUnitVector(i, 2)) for i in range(2)]
        msf2 = [basis2.makeFunction(ellipse, self.makeUnitVector(i, 3)) for i in range(3)]
        msf3 = [basis3.makeFunction(ellipse, self.makeUnitVector(i, 5)) for i in range(5)]
        for a, b in zip(msf3, msf1+msf2):
            self.compareMultiShapeletFunctions(a, b)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(MultiShapeletTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
