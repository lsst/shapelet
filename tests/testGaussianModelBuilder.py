#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

"""
Tests for GaussianModelBuilder

Run with:
   ./testGaussianModelBuilder.py
or
   python
   >>> import testGaussianModelBuilder; testGaussianModelBuilder.run()
"""

import unittest
import numpy

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions
import lsst.afw.geom as geom
import lsst.afw.geom.ellipses as ellipses
import lsst.afw.image
import lsst.afw.detection
import lsst.meas.extensions.multiShapelet as ms

numpy.random.seed(5)
numpy.set_printoptions(linewidth=120)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class GaussianModelBuilderTestCase(unittest.TestCase):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def buildModel(self, ellipse):
        model = numpy.zeros(self.region.getArea(), dtype=float).transpose()
        n = 0
        gt = ellipse.getGridTransform()
        for span in self.region.getSpans():
            y = span.getY()
            for x in range(span.getX0(), span.getX1() + 1):
                p = gt(geom.Point2D(x, y))
                model[n] = numpy.exp(-0.5 * (p.getX()**2 + p.getY()**2))
                n += 1
        return model

    def buildNumericalDerivative(self, builder, parameters, makeEllipse):
        eps = 1E-6
        derivative = numpy.zeros((len(parameters), builder.getSize()), dtype=float).transpose()
        for i in range(len(parameters)):
            parameters[i] += eps
            ellipse = makeEllipse(parameters)
            derivative[:,i] = builder.computeModel(ellipse)
            parameters[i] -= 2.0 * eps
            ellipse = makeEllipse(parameters)
            derivative[:,i] -= builder.computeModel(ellipse)
            derivative[:,i] /= 2.0 * eps
        return derivative

    def setUp(self):
        self.ellipse = ellipses.Ellipse(ellipses.Axes(10, 7, 0.3), geom.Point2D(500, 600))
        self.bbox = geom.Box2I(lsst.afw.geom.Point2I(480, 580), geom.Point2I(501, 601))
        self.region = lsst.afw.detection.Footprint(self.bbox)
        self.model = self.buildModel(self.ellipse)

    def tearDown(self):
        del self.ellipse
        del self.bbox
        del self.region

    def testModel1(self):
        builder = ms.GaussianModelBuilder(self.bbox)
        self.assertClose(builder.computeModel(self.ellipse), self.model)

    def testModel2(self):
        builder = ms.GaussianModelBuilder(self.region)
        self.assertClose(builder.computeModel(self.ellipse), self.model)

    def testDerivative1(self):
        """test derivative with no reparameterization"""
        builder = ms.GaussianModelBuilder(self.bbox)
        a = numpy.zeros((5, builder.getSize()), dtype=float).transpose()
        builder.computeDerivative(a, self.ellipse)
        def makeAxesEllipse(p):
            return ellipses.Ellipse(ellipses.Axes(*p[0:3]), geom.Point2D(*p[3:5]))
        n = self.buildNumericalDerivative(builder, self.ellipse.getParameterVector(), makeAxesEllipse)
        # no hard requirement for tolerances here, but I've dialed them to the max to avoid regressions
        self.assertClose(a, n, rtol=1E-4, atol=1E-6)

    def testDerivative2(self):
        """test derivative with trivial reparameterization (derivative wrt center point only)"""
        builder = ms.GaussianModelBuilder(self.bbox)
        jac = numpy.zeros((5, 2), dtype=float)
        jac[3,0] = 1.0
        jac[4,1] = 1.0
        a = numpy.zeros((2, builder.getSize()), dtype=float).transpose()
        builder.computeDerivative(a, self.ellipse, jac)
        def makePoint(p):
            return ellipses.Ellipse(self.ellipse.getCore(), geom.Point2D(*p))
        n = self.buildNumericalDerivative(builder, self.ellipse.getParameterVector()[3:5], makePoint)
        # no hard requirement for tolerances here, but I've dialed them to the max to avoid regressions
        self.assertClose(a, n, rtol=1E-4, atol=1E-6)

    def testDerivative3(self):
        """test derivative with nontrivial reparameterization (derivative wrt different core)"""
        builder = ms.GaussianModelBuilder(self.bbox)
        builder.computeModel(self.ellipse)
        quad = ellipses.Quadrupole(self.ellipse.getCore())
        jac = numpy.zeros((5, 3), dtype=float)
        jac[:3,:] = self.ellipse.getCore().dAssign(quad)
        a = numpy.zeros((3, builder.getSize()),dtype=float).transpose()
        builder.computeDerivative(a, self.ellipse, jac, False, True)
        def makeQuadrupole(p):
            return ellipses.Ellipse(ellipses.Quadrupole(*p), self.ellipse.getCenter())
        n = self.buildNumericalDerivative(builder, quad.getParameterVector(), makeQuadrupole)
        # no hard requirement for tolerances here, but I've dialed them to the max to avoid regressions
        self.assertClose(a, n, rtol=1E-4, atol=1E-6)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(GaussianModelBuilderTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
