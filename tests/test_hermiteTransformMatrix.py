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

try:
    import scipy.special
except ImportError:
    scipy = None

import lsst.utils.tests
import lsst.geom
import lsst.shapelet.tests


class HermiteTransformMatrixTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        np.random.seed(500)
        self.order = 4
        self.size = lsst.shapelet.computeSize(self.order)
        self.htm = lsst.shapelet.HermiteTransformMatrix(self.order)

    @staticmethod
    def ht(n):
        """return a scipy poly1d for the nth 'alternate' Hermite polynomial (i.e. Hermite polynomial
        with shapeley normalization)"""
        return (scipy.poly1d([(2**n * np.pi**0.5 * scipy.special.gamma(n+1))**(-0.5)])
                * scipy.special.hermite(n))

    def testCoefficientMatrices(self):
        coeff = self.htm.getCoefficientMatrix()
        coeffInv = self.htm.getInverseCoefficientMatrix()
        self.assertFloatsAlmostEqual(np.identity(self.order+1), np.dot(coeff, coeffInv))
        # Both matrices should be lower-triangular
        for i in range(0, self.order+1):
            for j in range(i+1, self.order+1):
                self.assertEqual(coeff[i, j], 0.0)
                self.assertEqual(coeffInv[i, j], 0.0)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def testCoefficientsAgainstHermite(self):
        """Test coefficient matrix values against scipy Hermite polynomials"""
        coeff = self.htm.getCoefficientMatrix()
        for n in range(0, self.order+1):
            poly = self.ht(n)
            self.assertFloatsAlmostEqual(coeff[n, :n+1], poly.c[::-1], atol=1E-15)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def testTransformMatrix(self):
        s = lsst.geom.LinearTransform.makeScaling(2.0, 1.5)
        r = lsst.geom.LinearTransform.makeRotation(0.30*lsst.geom.radians)
        transforms = [s, r, s*r*s]
        testPoints = np.random.randn(10, 2)
        for transform in transforms:
            m = self.htm.compute(transform)
            for testPoint in testPoints:
                assert testPoint.size == 2
                origPoint = lsst.geom.Point2D(testPoint[0], testPoint[1])
                transPoint = transform(origPoint)
                for i, inx, iny in lsst.shapelet.HermiteIndexGenerator(self.order):
                    v1 = self.ht(inx)(transPoint.getX()) * self.ht(iny)(transPoint.getY())
                    v2 = 0.0
                    for j, jnx, jny in lsst.shapelet.HermiteIndexGenerator(self.order):
                        v2 += m[i, j] * self.ht(jnx)(origPoint.getX()) * self.ht(jny)(origPoint.getY())
                    self.assertFloatsAlmostEqual(v1, v2, rtol=1E-11)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
