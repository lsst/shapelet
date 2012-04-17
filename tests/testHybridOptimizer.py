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
Tests for HybridOptimizer

Run with:
   ./testHybridOptimizer.py
or
   python
   >>> import testHybridOptimizer; testHybridOptimizer.run()
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
import testLib

numpy.random.seed(5)
numpy.set_printoptions(linewidth=120)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class HybridOptimizerTestCase(unittest.TestCase):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def testRosenbrock(self):
        ctrl = ms.HybridOptimizerControl()
        ctrl.fTol = 1E-10
        ctrl.gTol = 1E-10
        ctrl.minStep = 1E-14
        ctrl.maxIter = 200
        ctrl.tau = 1E-3
        initial = numpy.array([-1.2, 1.0], dtype=float)
        for lambda_ in (0.0, 1E-5, 1.0, 1E2, 1E4):
            obj = testLib.RosenbrockObjective(lambda_)
            opt = ms.HybridOptimizer(obj, initial, ctrl)
            for k in range(ctrl.maxIter):
                state = opt.step()
                if state & ms.HybridOptimizer.FINISHED:
                    break
            final = opt.getParameters()
            dx = ((final - numpy.array([1.0, 1.0], dtype=float))**2).sum()**0.5
            if lambda_ <= 1.0:
                self.assert_(state & ms.HybridOptimizer.SUCCESS)
            if state & ms.HybridOptimizer.SUCCESS_FTOL:
                self.assert_(numpy.abs(numpy.max(opt.getFunction())) <= ctrl.fTol)
            self.assert_(dx <= 1E-4)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(HybridOptimizerTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
