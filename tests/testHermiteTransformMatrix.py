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

try:
    import scipy.special
except ImportError:
    scipy = None

import lsst.utils.tests
import lsst.afw.geom
import lsst.shapelet.tests

numpy.random.seed(500)


class HermiteTransformMatrixTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        self.order = 4
        self.size = lsst.shapelet.computeSize(self.order)
        self.htm = lsst.shapelet.HermiteTransformMatrix(self.order)

    @staticmethod
    def ht(n):
        """return a scipy poly1d for the nth 'alternate' Hermite polynomial (i.e. Hermite polynomial
        with shapeley normalization)"""
        return (scipy.poly1d([(2**n * numpy.pi**0.5 * scipy.special.gamma(n+1))**(-0.5)])
                * scipy.special.hermite(n))

    def testCoefficientMatrices(self):
        coeff = self.htm.getCoefficientMatrix()
        coeffInv = self.htm.getInverseCoefficientMatrix()
        self.assertClose(numpy.identity(self.order+1), numpy.dot(coeff, coeffInv))
        # Both matrices should be lower-triangular
        for i in range(0, self.order+1):
            for j in range(i+1, self.order+1):
                self.assertEqual(coeff[i, j], 0.0)
                self.assertEqual(coeffInv[i, j], 0.0)
        # test coefficient matrix values against scipy Hermite polynomials
        if scipy is None:
            print "Skipping Hermite polynomial tests that require SciPy"
            return
        for n in range(0, self.order+1):
            poly = self.ht(n)
            self.assertClose(coeff[n, :n+1], poly.c[::-1], atol=1E-15)

    def testTransformMatrix(self):
        if scipy is None:
            print "Skipping transform tests that require SciPy"
            return

        s = lsst.afw.geom.LinearTransform.makeScaling(2.0, 1.5)
        r = lsst.afw.geom.LinearTransform.makeRotation(0.30*lsst.afw.geom.radians)
        transforms = [s, r, s*r*s]
        testPoints = numpy.random.randn(10, 2)
        for transform in transforms:
            m = self.htm.compute(transform)
            for testPoint in testPoints:
                assert(testPoint.size == 2)
                origPoint = lsst.afw.geom.Point2D(testPoint[0], testPoint[1])
                transPoint = transform(origPoint)
                for i, inx, iny in lsst.shapelet.HermiteIndexGenerator(self.order):
                    v1 = self.ht(inx)(transPoint.getX()) * self.ht(iny)(transPoint.getY())
                    v2 = 0.0
                    for j, jnx, jny in lsst.shapelet.HermiteIndexGenerator(self.order):
                        v2 += m[i, j] * self.ht(jnx)(origPoint.getX()) * self.ht(jny)(origPoint.getY())
                    self.assertClose(v1, v2, rtol=1E-11)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(HermiteTransformMatrixTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
