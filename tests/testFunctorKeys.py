from builtins import range
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

import lsst.utils.tests
import lsst.afw.geom.ellipses
import lsst.shapelet.tests
import lsst.afw.image
import lsst.afw.table

numpy.random.seed(500)


class FunctorKeyTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def testComputeOrder(self):
        invalidSizes = list(range(1, lsst.shapelet.computeSize(10)))
        for order in range(10):
            size = lsst.shapelet.computeSize(order)
            self.assertEqual(order, lsst.shapelet.computeOrder(size))
            invalidSizes.remove(size)
        for size in invalidSizes:
            self.assertRaises(lsst.pex.exceptions.InvalidParameterError, lsst.shapelet.computeOrder, size)

    def testShapeletFunctionKey(self):
        schema = lsst.afw.table.Schema()
        order = 4
        k0 = lsst.shapelet.ShapeletFunctionKey.addFields(schema, "s", "shapelet function",
                                                         "pixel", "count", order)
        k1 = lsst.shapelet.ShapeletFunctionKey(k0.getEllipse(), k0.getCoefficients())
        k2 = lsst.shapelet.ShapeletFunctionKey(schema["s"])
        self.assertEqual(k0, k1)
        self.assertEqual(k1, k2)
        self.assertTrue(k0.isValid())
        self.assertEqual(k0.getEllipse(), lsst.afw.table.EllipseKey(schema["s"]))
        self.assertEqual(k0.getCoefficients(), lsst.afw.table.ArrayDKey(schema["s"]))
        table = lsst.afw.table.BaseTable.make(schema)
        record = table.makeRecord()
        s0 = self.makeRandomShapeletFunction(order=order)
        record.set(k0, s0)
        s1 = record.get(k0)
        self.compareShapeletFunctions(s0, s1)
        self.assertRaises(lsst.pex.exceptions.InvalidParameterError, record.set, k0,
                          self.makeRandomShapeletFunction(order=3))

    def testMultiShapeletFunctionKey(self):
        schema = lsst.afw.table.Schema()
        msf0 = self.makeRandomMultiShapeletFunction(nComponents=3)
        orders = [s.getOrder() for s in msf0.getComponents()]
        k0 = lsst.shapelet.MultiShapeletFunctionKey.addFields(schema, "s", "shapelet function",
                                                              "pixel", "count", orders)
        k1 = lsst.shapelet.MultiShapeletFunctionKey([k0[i] for i in range(len(orders))])
        k2 = lsst.shapelet.MultiShapeletFunctionKey(schema["s"])
        self.assertEqual(k0, k1)
        self.assertEqual(k1, k2)
        self.assertTrue(k0.isValid())
        table = lsst.afw.table.BaseTable.make(schema)
        record = table.makeRecord()
        record.set(k0, msf0)
        msf1 = record.get(k0)
        self.compareMultiShapeletFunctions(msf0, msf1)
        self.assertRaises(lsst.pex.exceptions.InvalidParameterError, record.set, k0,
                          self.makeRandomMultiShapeletFunction(nComponents=4))
        self.assertRaises(lsst.pex.exceptions.NotFoundError,
                          lsst.shapelet.MultiShapeletFunctionKey,
                          schema["a"])


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()
    suites = []
    suites += unittest.makeSuite(FunctorKeyTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)


if __name__ == "__main__":
    run(True)
