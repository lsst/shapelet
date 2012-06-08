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
import lsst.afw.math

import scipy.ndimage
import numpy

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

    def checkConvolution(self, f1, f2):
        bbox = geom.Box2I(geom.Point2I(-50, -50), geom.Point2I(50, 50))
        i1 = lsst.afw.image.ImageD(bbox)
        f1.evaluate().addToImage(i1)
        self.assertClose(i1.getArray().sum(), f1.evaluate().integrate(), rtol=1E-3)
        i2 = lsst.afw.image.ImageD(bbox)
        f2.evaluate().addToImage(i2)
        self.assertClose(i2.getArray().sum(), f2.evaluate().integrate(), rtol=1E-3)
        fc1 = f1.convolve(f2)
        fc2 = f2.convolve(f1)
        ic1 = lsst.afw.image.ImageD(bbox)
        fc1.evaluate().addToImage(ic1)
        ic2 = lsst.afw.image.ImageD(bbox)
        fc2.evaluate().addToImage(ic2)
        self.assertClose(ic1.getArray(), ic2.getArray())
        out = lsst.afw.image.ImageD(bbox)
        # I'm using scipy.ndimage to convolve test images, because I can't figure
        # out how to make afw do it (afw can convolve images with kernels, but two similarly-sized
        # are apparently another matter; if I try to make a FixedKernel from one of the images,
        # I can't even make the operation commutative, let alone correct.
        scipy.ndimage.convolve(i1.getArray(), i2.getArray(), output=out.getArray(),
                               mode="constant", cval=0.0)
        self.assertClose(out.getArray(), ic1.getArray(), rtol=1E-4, atol=1E-5)
        self.assertClose(out.getArray(), ic2.getArray(), rtol=1E-4, atol=1E-5)
        return fc1, fc2

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
                p2 = numpy.dot(v, self.coefficients) * t.getLinear().computeDeterminant()
                self.assertClose(p1, p2)
            v = numpy.zeros(self.coefficients.shape, dtype=float)
            basis.fillIntegration(v)
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

    def testConvolution(self):
        e1 = ellipses.Ellipse(ellipses.Axes(10, 8, 0.3), geom.Point2D(1.5, 2.0))
        e2 = ellipses.Ellipse(ellipses.Axes(12, 9, -0.5), geom.Point2D(-1.0, -0.25))
        f1 = lsst.shapelet.ShapeletFunction(3, lsst.shapelet.HERMITE, e1)
        f2 = lsst.shapelet.ShapeletFunction(2, lsst.shapelet.LAGUERRE, e2)
        f1.getCoefficients()[:] = numpy.random.randn(*f1.getCoefficients().shape)
        f2.getCoefficients()[:] = numpy.random.randn(*f2.getCoefficients().shape)
        fc1, fc2 = self.checkConvolution(f1, f2)
        self.assertEqual(fc1.getBasisType(), lsst.shapelet.HERMITE)
        self.assertEqual(fc2.getBasisType(), lsst.shapelet.LAGUERRE)
        self.assertClose(fc1.getEllipse().getParameterVector(), fc2.getEllipse().getParameterVector())
        self.assertEqual(fc1.getOrder(), fc2.getOrder())
        fc2.changeBasisType(lsst.shapelet.HERMITE)
        self.assertClose(fc1.getCoefficients(), fc2.getCoefficients())
            
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
        model /= ellipse.getArea() / numpy.pi
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

    def testModel(self):
        builder = lsst.shapelet.ModelBuilder(self.x, self.y)
        builder.update(self.ellipse)
        z1 = numpy.random.randn(*self.model.transpose().shape).transpose()
        z0 = self.model + z1
        builder.addModelMatrix(self.order, z1)
        self.assertClose(z0, z1)
        coefficients = numpy.random.randn(self.model.shape[1])
        y1 = numpy.random.randn(self.model.shape[0])
        y0 = numpy.dot(self.model, coefficients) + y1
        builder.addModelVector(self.order, coefficients, y1)
        self.assertClose(y0, y1)

    def testMultiShapelet(self):
        """Should be redundant with testModel, but we want to be completely sure shapelet
        functions can be evaluated with ModelBuilder.addModelVector."""
        builder = lsst.shapelet.ModelBuilder(self.x, self.y)
        msf = lsst.shapelet.MultiShapeletFunction()
        z0 = numpy.zeros(self.model.shape[0], dtype=float)
        z1 = numpy.zeros(self.model.shape[0], dtype=float)
        for i in range(4):
            a, b = 6 * numpy.random.randn(2)**2
            theta = 3 * numpy.random.randn()
            axes = ellipses.Axes(a, b, theta)
            sf = lsst.shapelet.ShapeletFunction(
                2, lsst.shapelet.HERMITE,
                ellipses.Ellipse(axes, geom.Point2D(0,0)),
                numpy.random.randn(6)
                )
            msf.getElements().push_back(sf)
            builder.update(axes)
            builder.addModelVector(sf.getOrder(), sf.getCoefficients(), z1)
        ev = msf.evaluate()
        n = 0
        for x, y in zip(self.x, self.y):
            z0[n] = ev(x, y)
            n += 1
        self.assertClose(z1, z0)
            
            

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
