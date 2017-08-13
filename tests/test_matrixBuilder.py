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
from __future__ import absolute_import, division, print_function
from builtins import zip
import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.geom.ellipses
import lsst.shapelet.tests
import lsst.afw.image


class MatrixBuilderTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        np.random.seed(500)
        self.xD = np.random.randn(50)
        self.yD = np.random.randn(50)
        self.xF = self.xD.astype(np.float32)
        self.yF = self.yD.astype(np.float32)

    def checkAccessors(self, obj, basisSize):
        self.assertEqual(obj.getDataSize(), self.xD.size)
        self.assertEqual(obj.getBasisSize(), basisSize)

    def testShapeletMatrixBuilder(self):
        function = self.makeRandomShapeletFunction(order=4)
        size = function.getCoefficients().size
        function.getCoefficients()[:] = np.random.randn(size)
        basis = lsst.shapelet.MultiShapeletBasis(size)
        basis.addComponent(1.0, function.getOrder(), np.identity(size))
        factoryF = lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, function.getOrder())
        factoryD = lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, function.getOrder())
        # we should get the same results with an appropriately simple MultiShapeletBasis
        compoundFactoryF = lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, basis)
        compoundFactoryD = lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, basis)
        self.checkAccessors(factoryF, size)
        self.checkAccessors(factoryD, size)
        self.checkAccessors(compoundFactoryF, size)
        self.checkAccessors(compoundFactoryD, size)
        builder1F = factoryF()
        builder1D = factoryD()
        self.checkAccessors(builder1F, size)
        self.checkAccessors(builder1D, size)
        workspaceF = lsst.shapelet.MatrixBuilderF.Workspace(factoryF.computeWorkspace())
        workspaceD = lsst.shapelet.MatrixBuilderD.Workspace(factoryD.computeWorkspace())
        builder2F = factoryF(workspaceF)
        builder2D = factoryD(workspaceD)
        self.assertEqual(workspaceF.getRemaining(), 0)
        self.assertEqual(workspaceD.getRemaining(), 0)
        self.checkAccessors(builder2F, size)
        self.checkAccessors(builder2D, size)
        builder3F = compoundFactoryF()
        builder3D = compoundFactoryD()
        self.checkAccessors(builder3F, size)
        self.checkAccessors(builder3D, size)
        matrix1F = builder1F(function.getEllipse())
        matrix1D = builder1D(function.getEllipse())
        matrix2F = builder2F(function.getEllipse())
        matrix2D = builder2D(function.getEllipse())
        matrix3F = builder3F(function.getEllipse())
        matrix3D = builder3D(function.getEllipse())
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1D, matrix2D, rtol=0.0, atol=0.0)
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=0.0, atol=0.0)
        # same code, different construction pattern
        self.assertFloatsAlmostEqual(matrix1F, matrix3F, rtol=0.0, atol=0.0)
        # same code, different construction pattern
        self.assertFloatsAlmostEqual(matrix1D, matrix3D, rtol=0.0, atol=0.0)
        # same code, different construction pattern
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=1E-7, atol=0.0)
        # Finally, check against a completely different implementation (which is tested elsewhere)
        checkEvaluator = function.evaluate()
        checkVector = checkEvaluator(self.xD, self.yD)
        self.assertFloatsAlmostEqual(np.dot(matrix1D, function.getCoefficients()),
                                     checkVector, rtol=1E-13, atol=1E-14)

    def testRemappedShapeletMatrixBuilder(self):
        function = self.makeRandomShapeletFunction(order=4)
        size = 6
        radius = 3.2
        remapMatrix = np.random.randn(function.getCoefficients().size, size)
        coefficients = np.random.randn(size)
        function.getCoefficients()[:] = np.dot(remapMatrix, coefficients)
        basis = lsst.shapelet.MultiShapeletBasis(size)
        basis.addComponent(radius, function.getOrder(), remapMatrix)
        factoryF = lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, basis)
        factoryD = lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, basis)
        self.checkAccessors(factoryF, size)
        self.checkAccessors(factoryD, size)
        builder1F = factoryF()
        builder1D = factoryD()
        self.checkAccessors(builder1F, size)
        self.checkAccessors(builder1D, size)
        workspaceF = lsst.shapelet.MatrixBuilderF.Workspace(factoryF.computeWorkspace())
        workspaceD = lsst.shapelet.MatrixBuilderD.Workspace(factoryD.computeWorkspace())
        builder2F = factoryF(workspaceF)
        builder2D = factoryD(workspaceD)
        self.assertEqual(workspaceF.getRemaining(), 0)
        self.assertEqual(workspaceD.getRemaining(), 0)
        self.checkAccessors(builder2F, size)
        self.checkAccessors(builder2D, size)
        matrix1F = builder1F(function.getEllipse())
        matrix1D = builder1D(function.getEllipse())
        matrix2F = builder2F(function.getEllipse())
        matrix2D = builder2D(function.getEllipse())
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1D, matrix2D, rtol=0.0, atol=0.0)
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=0.0, atol=0.0)
        # same code, different precision
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=1E-7, atol=0.0)
        # Finally, check against a completely different implementation (which is tested elsewhere)
        function.getEllipse().scale(radius)
        checkEvaluator = function.evaluate()
        checkVector = checkEvaluator(self.xD, self.yD)
        self.assertFloatsAlmostEqual(np.dot(matrix1D, coefficients), checkVector, rtol=1E-13, atol=1E-14)

    def testConvolvedShapeletMatrixBuilder(self):
        function = self.makeRandomShapeletFunction(order=4)
        psf = self.makeRandomMultiShapeletFunction(nComponents=1)
        size = function.getCoefficients().size
        function.getCoefficients()[:] = np.random.randn(size)
        basis = lsst.shapelet.MultiShapeletBasis(size)
        basis.addComponent(1.0, function.getOrder(), np.identity(size))
        factoriesF = [lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, function.getOrder(),
                                                           psf.getComponents()[0]),
                      lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, basis, psf),
                      ]
        factoriesD = [lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, function.getOrder(),
                                                           psf.getComponents()[0]),
                      lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, basis, psf),
                      ]
        lastF = None
        for factory in factoriesF:
            self.checkAccessors(factory, size)
            workspace = lsst.shapelet.MatrixBuilderF.Workspace(factory.computeWorkspace())
            builder1 = factory()
            builder2 = factory(workspace)
            self.checkAccessors(builder1, size)
            self.checkAccessors(builder2, size)
            self.assertEqual(workspace.getRemaining(), 0)
            matrix1 = builder1(function.getEllipse())
            matrix2 = builder2(function.getEllipse())
            # same code, different workspace
            self.assertFloatsAlmostEqual(matrix1, matrix2, rtol=0.0, atol=0.0)
            if lastF is not None:
                # same code, different construction
                self.assertFloatsAlmostEqual(matrix1, lastF, rtol=0.0, atol=0.0)
            lastF = matrix1
        lastD = None
        for factory in factoriesD:
            self.checkAccessors(factory, size)
            workspace = lsst.shapelet.MatrixBuilderD.Workspace(factory.computeWorkspace())
            builder1 = factory()
            builder2 = factory(workspace)
            self.checkAccessors(builder1, size)
            self.checkAccessors(builder2, size)
            self.assertEqual(workspace.getRemaining(), 0)
            matrix1 = builder1(function.getEllipse())
            matrix2 = builder2(function.getEllipse())
            # same code, different workspace
            self.assertFloatsAlmostEqual(matrix1, matrix2, rtol=0.0, atol=0.0)
            if lastD is not None:
                # same code, different construction
                self.assertFloatsAlmostEqual(matrix1, lastD, rtol=0.0, atol=0.0)
            lastD = matrix1
        # same code, different precision
        self.assertFloatsAlmostEqual(lastF, lastD, rtol=1E-6, atol=1E-5)

        # Finally, check against a completely different implementation (which is tested elsewhere)
        convolved = function.convolve(psf.getComponents()[0])
        checkEvaluator = convolved.evaluate()
        checkVector = checkEvaluator(self.xD, self.yD)
        self.assertFloatsAlmostEqual(np.dot(lastD, function.getCoefficients()), checkVector, rtol=1E-13)

    def testRemappedConvolvedShapeletMatrixBuilder(self):
        function = self.makeRandomShapeletFunction(order=4)
        psf = self.makeRandomMultiShapeletFunction(nComponents=1)
        size = 6
        radius = 3.2
        remapMatrix = np.random.randn(function.getCoefficients().size, size)
        coefficients = np.random.randn(size)
        function.getCoefficients()[:] = np.dot(remapMatrix, coefficients)
        basis = lsst.shapelet.MultiShapeletBasis(size)
        basis.addComponent(radius, function.getOrder(), remapMatrix)
        factoryF = lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, basis, psf)
        factoryD = lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, basis, psf)
        matrixF = None
        for factory in (factoryF, factoryD):
            self.checkAccessors(factory, size)
            workspace = factory.Workspace(factory.computeWorkspace())
            builder1 = factory()
            builder2 = factory(workspace)
            self.checkAccessors(builder1, size)
            self.checkAccessors(builder2, size)
            self.assertEqual(workspace.getRemaining(), 0)
            matrix1 = builder1(function.getEllipse())
            matrix2 = builder2(function.getEllipse())
            # same code, different workspace
            self.assertFloatsAlmostEqual(matrix1, matrix2, rtol=0.0, atol=0.0)
            if matrixF is None:
                matrixF = matrix1
        matrixD = matrix1
        # same code, different precision
        self.assertFloatsAlmostEqual(matrixF, matrixD, rtol=1E-4, atol=0.0)

        # Finally, check against a completely different implementation (which is tested elsewhere)
        function.getEllipse().scale(radius)
        convolved = function.convolve(psf.getComponents()[0])
        checkEvaluator = convolved.evaluate()
        checkVector = checkEvaluator(self.xD, self.yD)
        self.assertFloatsAlmostEqual(np.dot(matrixD, coefficients), checkVector, rtol=1E-12, atol=1E-12)

    def testCompoundMatrixBuilder(self):
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(4.0, 3.0, 1.0),
                                                 lsst.afw.geom.Point2D(3.2, 1.0))
        radii = [0.7, 1.2]
        orders = [4, 3]
        size = 8
        functions = [self.makeRandomShapeletFunction(order=order, ellipse=ellipse, scale=radius)
                     for radius, order in zip(radii, orders)]
        remapMatrices = [np.random.randn(f.getCoefficients().size, size) for f in functions]
        coefficients = np.random.randn(size)
        basis = lsst.shapelet.MultiShapeletBasis(size)
        for radius, function, matrix in zip(radii, functions, remapMatrices):
            function.getCoefficients()[:] = np.dot(matrix, coefficients)
            basis.addComponent(radius, function.getOrder(), matrix)
        factoryF = lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, basis)
        factoryD = lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, basis)
        self.checkAccessors(factoryF, size)
        self.checkAccessors(factoryD, size)
        builder1F = factoryF()
        builder1D = factoryD()
        self.checkAccessors(builder1F, size)
        self.checkAccessors(builder1D, size)
        workspaceF = lsst.shapelet.MatrixBuilderF.Workspace(factoryF.computeWorkspace())
        workspaceD = lsst.shapelet.MatrixBuilderD.Workspace(factoryD.computeWorkspace())
        builder2F = factoryF(workspaceF)
        builder2D = factoryD(workspaceD)
        self.assertEqual(workspaceF.getRemaining(), 0)
        self.assertEqual(workspaceD.getRemaining(), 0)
        self.checkAccessors(builder2F, size)
        self.checkAccessors(builder2D, size)
        matrix1F = builder1F(ellipse)
        matrix1D = builder1D(ellipse)
        matrix2F = builder2F(ellipse)
        matrix2D = builder2D(ellipse)
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1D, matrix2D, rtol=0.0, atol=0.0)
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=0.0, atol=0.0)
        # same code, different precision
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=1E-7, atol=0.0)
        # Finally, check against a completely different implementation (which is tested elsewhere)
        checkVector = np.zeros(self.xD.shape, dtype=float)
        for function in functions:
            checkEvaluator = function.evaluate()
            checkVector += checkEvaluator(self.xD, self.yD)
        self.assertFloatsAlmostEqual(np.dot(matrix1D, coefficients), checkVector, rtol=1E-13, atol=1E-14)

    def testConvolvedCompoundMatrixBuilder(self):
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(4.0, 3.0, 1.0),
                                                 lsst.afw.geom.Point2D(3.2, 1.0))
        radii = [0.7, 1.2]
        orders = [4, 3]
        size = 8
        functions = [self.makeRandomShapeletFunction(order=order, ellipse=ellipse, scale=radius)
                     for radius, order in zip(radii, orders)]
        psf = self.makeRandomMultiShapeletFunction()
        remapMatrices = [np.random.randn(f.getCoefficients().size, size) for f in functions]
        coefficients = np.random.randn(size)
        basis = lsst.shapelet.MultiShapeletBasis(size)
        for radius, function, matrix in zip(radii, functions, remapMatrices):
            function.getCoefficients()[:] = np.dot(matrix, coefficients)
            basis.addComponent(radius, function.getOrder(), matrix)
        msf = lsst.shapelet.MultiShapeletFunction(functions)
        factoryF = lsst.shapelet.MatrixBuilderF.Factory(self.xF, self.yF, basis, psf)
        factoryD = lsst.shapelet.MatrixBuilderD.Factory(self.xD, self.yD, basis, psf)
        self.checkAccessors(factoryF, size)
        self.checkAccessors(factoryD, size)
        builder1F = factoryF()
        builder1D = factoryD()
        self.checkAccessors(builder1F, size)
        self.checkAccessors(builder1D, size)
        workspaceF = lsst.shapelet.MatrixBuilderF.Workspace(factoryF.computeWorkspace())
        workspaceD = lsst.shapelet.MatrixBuilderD.Workspace(factoryD.computeWorkspace())
        builder2F = factoryF(workspaceF)
        builder2D = factoryD(workspaceD)
        self.assertEqual(workspaceF.getRemaining(), 0)
        self.assertEqual(workspaceD.getRemaining(), 0)
        self.checkAccessors(builder2F, size)
        self.checkAccessors(builder2D, size)
        matrix1F = builder1F(ellipse)
        matrix1D = builder1D(ellipse)
        matrix2F = builder2F(ellipse)
        matrix2D = builder2D(ellipse)
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1D, matrix2D, rtol=0.0, atol=0.0)
        # same code, different workspace
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=0.0, atol=0.0)
        # same code, different precision
        self.assertFloatsAlmostEqual(matrix1F, matrix2F, rtol=1E-7, atol=0.0)

        # Finally, check against a completely different implementation (which is tested elsewhere)
        convolved = msf.convolve(psf)
        checkEvaluator = convolved.evaluate()
        checkVector = checkEvaluator(self.xD, self.yD)
        self.assertFloatsAlmostEqual(np.dot(matrix1D, coefficients), checkVector, rtol=1E-13, atol=1E-14)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
