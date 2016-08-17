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
import cPickle

try:
    import scipy.special
except ImportError:
    scipy = None

import lsst.utils.tests
import lsst.shapelet.tests
import lsst.afw.image
import lsst.afw.geom as geom
import lsst.afw.geom.ellipses as ellipses

numpy.random.seed(500)


class ShapeletFunctionTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        order = 4
        self.ellipse = ellipses.Ellipse(ellipses.Axes(2.2, 0.8, 0.3), geom.Point2D(0.12, -0.08))
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

    def testPickle(self):
        for function in self.functions:
            s = cPickle.dumps(function, protocol=2)
            function2 = cPickle.loads(s)
            self.assertEqual(function.getOrder(), function2.getOrder())
            self.assertEqual(function.getBasisType(), function2.getBasisType())
            self.assertClose(function.getEllipse().getParameterVector(),
                             function2.getEllipse().getParameterVector())
            self.assertClose(function.getCoefficients(), function2.getCoefficients())

    def testConversion(self):
        for basis, function in zip(self.bases, self.functions):
            evaluator = function.evaluate()
            v = numpy.zeros(self.coefficients.shape, dtype=float)
            t = self.ellipse.getGridTransform()
            for x, y in zip(self.x, self.y):
                basis.fillEvaluation(v, t(geom.Point2D(x, y)))
                p1 = evaluator(x, y)
                p2 = numpy.dot(v, self.coefficients) * t.getLinear().computeDeterminant()
                self.assertClose(p1, p2, rtol=1E-8)
            v = numpy.zeros(self.coefficients.shape, dtype=float)
            basis.fillIntegration(v)
            p1 = evaluator.integrate()
            p2 = numpy.dot(v, self.coefficients)
            self.assertClose(p1, p2, rtol=1E-8)

    def testMoments(self):
        x = numpy.linspace(-15, 15, 151)
        y = x
        for function in self.functions:
            z = self.makeImage(function, x, y)
            self.checkMoments(function, x, y, z)

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
                self.assertClose(dx_n, dx_a, rtol=1E-5)
                self.assertClose(dy_n, dy_a, rtol=1E-5)

    def testAddToImage(self):
        bbox = geom.Box2I(geom.Point2I(5, 6), geom.Extent2I(20, 30))
        image = lsst.afw.image.ImageD(bbox)
        x = numpy.arange(bbox.getBeginX(), bbox.getEndX(), dtype=float)
        y = numpy.arange(bbox.getBeginY(), bbox.getEndY(), dtype=float)
        array = numpy.zeros((bbox.getHeight(), bbox.getWidth()), dtype=float)
        for f in self.functions:
            image.getArray()[:] = 0.0
            array[:] = 0.0
            ev = f.evaluate()
            ev.addToImage(image)
            ev.addToImage(array, bbox.getMin())
            check = self.makeImage(f, x, y)
            self.assertClose(image.getArray(), check)
            self.assertClose(array, check)

    def testConvolution(self):
        if scipy is None:
            print "Skipping convolution test; scipy could not be imported."
            return
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
        self.assertClose(fc1.getCoefficients(), fc2.getCoefficients(), 1E-8)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ShapeletFunctionTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
