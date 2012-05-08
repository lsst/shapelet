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
Tests for math.shapelets

Run with:
   ./shapelets.py
or
   python
   >>> import shapelets; shapelets.run()
"""

import unittest
import numpy

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions
import lsst.afw.geom as geom
import lsst.afw.geom.ellipses as ellipses
import lsst.shapelet
import lsst.afw.image
import lsst.afw.detection

numpy.random.seed(5)
numpy.set_printoptions(linewidth=120)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ShapeletTestMixin(object):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def makeImage(self, function, x, y):
        z = numpy.zeros((y.size, x.size), dtype=float)
        e = function.evaluate()
        for i, py in enumerate(y):
            for j, px in enumerate(x):
                z[i,j] = e(float(px), float(py))
        return z

    def measureMoments(self, function, x, y, z):
        gx, gy = numpy.meshgrid(x, y)
        m = z.sum()
        dipole = geom.Point2D((gx * z).sum() / m, (gy * z).sum() / m)
        gx -= dipole.getX()
        gy -= dipole.getY()
        quadrupole = ellipses.Quadrupole(
            (gx**2 * z).sum() / m,
            (gy**2 * z).sum() / m,
            (gx * gy * z).sum() / m
            )
        imageMoments = ellipses.Ellipse(quadrupole, dipole)
        shapeletMoments = function.evaluate().computeMoments()
        self.assertClose(imageMoments.getCenter().getX(), shapeletMoments.getCenter().getX(),
                         rtol=1E-4, atol=1E-3)
        self.assertClose(imageMoments.getCenter().getY(), shapeletMoments.getCenter().getY(),
                         rtol=1E-4, atol=1E-3)
        self.assertClose(imageMoments.getCore().getIxx(), shapeletMoments.getCore().getIxx(),
                         rtol=1E-4, atol=1E-3)
        self.assertClose(imageMoments.getCore().getIyy(), shapeletMoments.getCore().getIyy(),
                         rtol=1E-4, atol=1E-3)
        self.assertClose(imageMoments.getCore().getIxy(), shapeletMoments.getCore().getIxy(),
                         rtol=1E-4, atol=1E-3)
        integral = numpy.trapz(numpy.trapz(z, gx, axis=1), y, axis=0)
        self.assertClose(integral, function.evaluate().integrate(), rtol=1E-3, atol=1E-2)

class ShapeletTestCase(unittest.TestCase, ShapeletTestMixin):

    def setUp(self):
        order = 4
        self.ellipse = ellipses.Ellipse(ellipses.Axes(1.2, 0.8, 0.3), geom.Point2D(0.12, -0.08))
        self.coefficients = numpy.random.randn(lsst.shapelet.computeSize(order))
        self.x = numpy.random.randn(25)
        self.y = numpy.random.randn(25)
        self.bases = [
            lsst.shapelet.BasisEvaluator(order, lsst.shapelet.HERMITE),
            lsst.shapelet.BasisEvaluator(order, lsst.shapelet.LAGUERRE),
            ]
        self.functions = [
            lsst.shapelet.ShapeletFunction(order, lsst.shapelet.HERMITE, self.coefficients),
            lsst.shapelet.ShapeletFunction(order, lsst.shapelet.LAGUERRE, self.coefficients),
            ]
        for function in self.functions:
            function.setEllipse(self.ellipse)

    def testConversion(self):
        for basis, function in zip(self.bases, self.functions):
            evaluator = function.evaluate()
            v = numpy.zeros(self.coefficients.shape, dtype=float)
            t = self.ellipse.getGridTransform()
            for x, y in zip(self.x, self.y):
                basis.fillEvaluation(v, t(geom.Point2D(x, y)))
                p1 = evaluator(x, y)
                p2 = numpy.dot(v, self.coefficients)
                self.assertClose(p1, p2)
            v = numpy.zeros(self.coefficients.shape, dtype=float)
            basis.fillIntegration(v)
            v /= t.getLinear().computeDeterminant()
            p1 = evaluator.integrate()
            p2 = numpy.dot(v, self.coefficients)
            self.assertClose(p1, p2)

    def testMoments(self):
        x = numpy.linspace(-15, 15, 151)
        y = x
        for function in self.functions:
            z = self.makeImage(function, x, y)
            self.measureMoments(function, x, y, z)

    def testDerivatives(self):
        eps = 1E-7
        v = numpy.zeros(self.coefficients.shape, dtype=float)
        v_lo = numpy.zeros(self.coefficients.shape, dtype=float)
        v_hi = numpy.zeros(self.coefficients.shape, dtype=float)
        dx_a = numpy.zeros(self.coefficients.shape, dtype=float)
        dy_a = numpy.zeros(self.coefficients.shape, dtype=float)
        for basis in self.bases:
            for x, y in zip(self.x, self.y):
                basis.fillEvaluation(v, x, y, dx_a, dy_a)
                basis.fillEvaluation(v_hi, x+eps, y) 
                basis.fillEvaluation(v_lo, x-eps, y)
                dx_n = 0.5 * (v_hi - v_lo) / eps
                basis.fillEvaluation(v_hi, x, y+eps) 
                basis.fillEvaluation(v_lo, x, y-eps)
                dy_n = 0.5 * (v_hi - v_lo) / eps
                self.assertClose(dx_n, dx_a, rtol=2.0*eps)
                self.assertClose(dy_n, dy_a, rtol=2.0*eps)

    def testAddToImage(self):
        bbox = geom.Box2I(geom.Point2I(5, 6), geom.Extent2I(20, 30))
        image = lsst.afw.image.ImageD(bbox)
        x = numpy.arange(bbox.getBeginX(), bbox.getEndX())
        y = numpy.arange(bbox.getBeginY(), bbox.getEndY())
        array = numpy.zeros((bbox.getHeight(), bbox.getWidth()), dtype=float)
        for f in self.functions:
            ev = f.evaluate()
            ev.addToImage(image)
            ev.addToImage(array, bbox.getMin())
            check = self.makeImage(f, x, y)
            self.assertClose(image.getArray(), check)
            self.assertClose(array, check)

class MultiShapeletTestCase(unittest.TestCase, ShapeletTestMixin):

    def testMoments(self):
        x = numpy.linspace(-50, 50, 1001)
        y = x
        elements = []
        for n in range(3):
            ellipse = ellipses.Ellipse(
                ellipses.Axes(
                    float(numpy.random.uniform(low=1, high=2)),
                    float(numpy.random.uniform(low=1, high=2)),
                    float(numpy.random.uniform(low=0, high=numpy.pi))
                    ),
                geom.Point2D(0.23, -0.15)
                )
            coefficients = numpy.random.randn(lsst.shapelet.computeSize(n))
            element = lsst.shapelet.ShapeletFunction(n, lsst.shapelet.HERMITE, coefficients)
            element.setEllipse(ellipse)
            elements.append(element)
        function = lsst.shapelet.MultiShapeletFunction(elements)
        x = numpy.linspace(-10, 10, 101)
        y = x
        z = self.makeImage(function, x, y)
        self.measureMoments(function, x, y, z)

class ModelBuilderTestCase(unittest.TestCase, ShapeletTestMixin):

    def buildModel(self, ellipse):
        model = numpy.zeros((lsst.shapelet.computeSize(self.order), self.x.size), dtype=float).transpose()
        evaluator = lsst.shapelet.BasisEvaluator(self.order, lsst.shapelet.HERMITE)
        n = 0
        gt = ellipse.getGridTransform()
        for x, y in zip(self.x, self.y):
            p = gt(geom.Point2D(x, y))
            evaluator.fillEvaluation(model[n,:], p)
            n += 1
        return model

    def setUp(self):
        self.order = 3
        self.ellipse = ellipses.Axes(10, 7, 0.3)
        self.xg, self.yg = numpy.meshgrid(numpy.linspace(-20, 20, 101), numpy.linspace(-15, 25, 95))
        self.x = self.xg.ravel()
        self.y = self.yg.ravel()
        self.model = self.buildModel(self.ellipse)

    def tearDown(self):
        del self.ellipse

    def check(self, builder, region, model):
        self.assertClose(builder.getModel(), model)

    def testModel(self):
        builder = lsst.shapelet.ModelBuilder(self.order, self.x, self.y)
        builder.update(self.ellipse)
        self.assertClose(builder.getModel(), self.model)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ShapeletTestCase)
    suites += unittest.makeSuite(MultiShapeletTestCase)
    suites += unittest.makeSuite(ModelBuilderTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
