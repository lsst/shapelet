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

import lsst.utils.tests
import lsst.afw.geom.ellipses
import lsst.shapelet.tests
import lsst.afw.image

numpy.random.seed(500)

class ModelBuilderTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def buildModel(self, ellipse):
        model = numpy.zeros((lsst.shapelet.computeSize(self.order), self.x.size), dtype=float).transpose()
        evaluator = lsst.shapelet.BasisEvaluator(self.order, lsst.shapelet.HERMITE)
        n = 0
        gt = ellipse.getGridTransform()
        for x, y in zip(self.x, self.y):
            p = gt(lsst.afw.geom.Point2D(x, y))
            evaluator.fillEvaluation(model[n,:], p)
            n += 1
        model /= ellipse.getArea() / numpy.pi
        return model

    def setUp(self):
        self.order = 3
        self.ellipse = lsst.afw.geom.ellipses.Axes(10, 7, 0.3)
        self.xg, self.yg = numpy.meshgrid(numpy.linspace(-20, 20, 101), numpy.linspace(-15, 25, 95))
        self.x = self.xg.ravel()
        self.y = self.yg.ravel()
        self.model = self.buildModel(self.ellipse)

    def tearDown(self):
        del self.ellipse

    def testModel(self):
        builder = lsst.shapelet.ModelBuilderD(self.x, self.y)
        builder.update(self.ellipse)
        z1 = numpy.random.randn(*self.model.transpose().shape).transpose()
        z0 = self.model + z1
        builder.addModelMatrix(self.order, z1)
        self.assertClose(z0, z1)
        coefficients = numpy.random.randn(self.model.shape[1])
        y1 = numpy.random.randn(self.model.shape[0])
        y0 = numpy.dot(self.model, coefficients) + y1
        builder.addModelVector(self.order, coefficients, y1)
        self.assertClose(y0, y1, rtol=1E-12)

    def testMultiShapeletFunction(self):
        """Should be redundant with testModel, but we want to be completely sure shapelet
        functions can be evaluated with ModelBuilder.addModelVector."""
        builder = lsst.shapelet.ModelBuilderD(self.x, self.y)
        msf = lsst.shapelet.MultiShapeletFunction()
        z0 = numpy.zeros(self.model.shape[0], dtype=float)
        z1 = numpy.zeros(self.model.shape[0], dtype=float)
        for i in range(4):
            a, b = 6 * numpy.random.randn(2)**2
            theta = 3 * numpy.random.randn()
            axes = lsst.afw.geom.ellipses.Axes(a, b, theta)
            sf = lsst.shapelet.ShapeletFunction(
                2, lsst.shapelet.HERMITE,
                lsst.afw.geom.ellipses.Ellipse(axes, lsst.afw.geom.Point2D(0,0)),
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
        self.assertClose(z1, z0, rtol=1E-8)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ModelBuilderTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
