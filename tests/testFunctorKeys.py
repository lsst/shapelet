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

numpy.random.seed(500)

class FunctorKeyTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def testComputeOrder(self):
        invalidSizes = range(1, lsst.shapelet.computeSize(10))
        for order in range(10):
            size = lsst.shapelet.computeSize(order)
            self.assertEqual(order, lsst.shapelet.computeOrder(size))
            invalidSizes.remove(size)
        for size in invalidSizes:
            self.assertRaises(lsst.pex.exceptions.InvalidParameterError, lsst.shapelet.computeOrder, size)

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
