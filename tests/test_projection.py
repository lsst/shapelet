#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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

import numpy as np

import lsst.geom
import lsst.utils.tests
import lsst.pex.exceptions
import lsst.afw.geom.ellipses
import lsst.shapelet.tests
import lsst.afw.image


class ProjectionTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        np.random.seed(500)
        self.ghp = lsst.shapelet.GaussHermiteProjection(16)

    def testRotation(self):
        order = 4
        size = lsst.shapelet.computeSize(order)
        nPoints = 100
        unitCircle = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Quadrupole(1.0, 1.0, 0.0))
        # This matrix should represent a pure rotation in shapelet space, which can be done exactly.
        inputTransform = lsst.geom.LinearTransform.makeRotation(0.0*lsst.geom.radians)
        outputTransform = lsst.geom.LinearTransform.makeRotation(np.pi/3*lsst.geom.radians)
        m = self.ghp.compute(inputTransform, order, outputTransform, order)
        # If we apply a rotation by np.pi/3 six times, we should get back where we started with.
        self.assertFloatsAlmostEqual(np.linalg.matrix_power(m, 6),
                                     np.identity(size), rtol=1E-14, atol=1E-14)
        # Now we test that we get the same result (up to round-off error) for a bunch of test points.
        inputTestPoints = np.random.randn(2, nPoints)
        outputTestPoints = np.dot(outputTransform.getMatrix(), inputTestPoints)
        inputBuilder = lsst.shapelet.MatrixBuilderD.Factory(inputTestPoints[0, :],
                                                            inputTestPoints[1, :], order)()
        outputBuilder = lsst.shapelet.MatrixBuilderD.Factory(outputTestPoints[0, :],
                                                             outputTestPoints[1, :], order)()
        inputBasis = inputBuilder(unitCircle)
        outputBasis = outputBuilder(unitCircle)
        self.assertFloatsAlmostEqual(np.dot(outputBasis, m), inputBasis, rtol=1E-13)

    def testPerturbation(self):
        inputEllipse = lsst.afw.geom.ellipses.Quadrupole(2.0, 2.0, 0.0)
        outputEllipse = lsst.afw.geom.ellipses.Quadrupole(2.01, 2.0, 0.02)
        inputShapelet = lsst.shapelet.ShapeletFunction(2, lsst.shapelet.HERMITE,
                                                       lsst.afw.geom.ellipses.Ellipse(inputEllipse))
        outputShapelet = lsst.shapelet.ShapeletFunction(8, lsst.shapelet.HERMITE,
                                                        lsst.afw.geom.ellipses.Ellipse(outputEllipse))
        m = self.ghp.compute(inputEllipse, inputShapelet.getOrder(),
                             outputEllipse, outputShapelet.getOrder())
        x = np.random.randn(50)
        y = np.random.randn(50)
        for i, nx, ny in lsst.shapelet.HermiteIndexGenerator(2):
            inputShapelet.getCoefficients()[:] = 0.0
            inputShapelet.getCoefficients()[i] = 1.0
            outputShapelet.getCoefficients()[:] = np.dot(m, inputShapelet.getCoefficients())
            inputEv = inputShapelet.evaluate()
            outputEv = outputShapelet.evaluate()
            inputZ = np.array([inputEv(px, py) for px, py in zip(x, y)])
            outputZ = np.array([outputEv(px, py) for px, py in zip(x, y)])
            # tolerances here aren't rigorous, because we're testing an approximation without
            # computing how good it should be; it's mostly just a regression/sanity test
            self.assertFloatsAlmostEqual(inputZ, outputZ, rtol=5E-3)

    def testConvolution0(self):
        """Test that the specialization for zeroth-order convolution produces the same tensor
        as the general case, which is tested more rigorously elsewhere."""
        psf = self.makeRandomShapeletFunction(order=4)
        ellipse1 = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Quadrupole(2.0, 3.0, 1.0))
        ellipse2 = lsst.afw.geom.ellipses.Ellipse(ellipse1)
        ghcN = lsst.shapelet.GaussHermiteConvolution(1, psf)
        ghc0 = lsst.shapelet.GaussHermiteConvolution(0, psf)
        rN = ghcN.evaluate(ellipse1)
        r0 = ghc0.evaluate(ellipse2)
        self.assertFloatsAlmostEqual(rN[:15, 0], r0[:15, 0], rtol=1E-14)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
