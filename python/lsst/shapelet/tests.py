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

"""
Test utility code for shapelets library; here so it can be used
in multiple test scripts and tests in downstream packages.
"""

import numpy
try:
    import scipy.ndimage
except ImportError:
    scipy = None

import lsst.utils.tests
import lsst.afw.geom
import lsst.afw.geom.ellipses


class ShapeletTestCase(lsst.utils.tests.TestCase):

    @staticmethod
    def makeUnitVector(i, n):
        v = numpy.zeros(n, dtype=float)
        v[i] = 1.0
        return v

    @staticmethod
    def makeImage(function, x, y):
        z = numpy.zeros((y.size, x.size), dtype=float)
        e = function.evaluate()
        for i, py in enumerate(y):
            for j, px in enumerate(x):
                z[i, j] = e(float(px), float(py))
        return z

    @staticmethod
    def makeRandomShapeletFunction(order=2, zeroCenter=False, ellipse=None, scale=1.0):
        center = lsst.afw.geom.Point2D()
        if not zeroCenter:
            center = lsst.afw.geom.Point2D(numpy.random.randn(), numpy.random.randn())
        if ellipse is None:
            ellipse = lsst.afw.geom.ellipses.Ellipse(
                lsst.afw.geom.ellipses.Axes(
                    float(numpy.random.uniform(low=1, high=2)),
                    float(numpy.random.uniform(low=1, high=2)),
                    float(numpy.random.uniform(low=0, high=numpy.pi))
                ),
                center
            )
        coefficients = numpy.random.randn(lsst.shapelet.computeSize(order))
        result = lsst.shapelet.ShapeletFunction(order, lsst.shapelet.HERMITE, coefficients)
        result.setEllipse(ellipse)
        result.getEllipse().scale(scale)
        return result

    @staticmethod
    def makeRandomMultiShapeletFunction(nComponents=3, ellipse=None):
        components = []
        for n in range(nComponents):
            components.append(ShapeletTestCase.makeRandomShapeletFunction(ellipse=ellipse))
        return lsst.shapelet.MultiShapeletFunction(components)

    def compareShapeletFunctions(self, a, b, rtolEllipse=1E-13, rtolCoeff=1E-13,
                                 atolEllipse=1E-14, atolCoeff=1E-14):
        self.assertEqual(a.getOrder(), b.getOrder())
        self.assertEqual(a.getBasisType(), b.getBasisType())
        self.assertClose(a.getEllipse().getParameterVector(), b.getEllipse().getParameterVector(),
                         rtol=rtolEllipse, atol=atolEllipse)
        self.assertClose(a.getCoefficients(), b.getCoefficients(), rtol=rtolCoeff, atol=atolCoeff)

    def simplifyMultiShapeletFunction(self, msf):
        keep = []
        for s in msf.getComponents():
            if not numpy.allclose(s.getCoefficients(), 0.0):
                params = tuple(s.getEllipse().getParameterVector()) + tuple(s.getCoefficients())
                keep.append((params, s))
        msf = lsst.shapelet.MultiShapeletFunction()
        keep.sort(key=lambda t: t[0])
        for params, s in keep:
            msf.getComponents().push_back(s)
        return msf

    def compareMultiShapeletFunctions(self, a, b, simplify=True, rtolEllipse=1E-13, rtolCoeff=1E-13,
                                      atolEllipse=1E-14, atolCoeff=1E-14):
        if simplify:
            a = self.simplifyMultiShapeletFunction(a)
            b = self.simplifyMultiShapeletFunction(b)
        self.assertEqual(a.getComponents().size(), b.getComponents().size())
        for sa, sb in zip(a.getComponents(), b.getComponents()):
            self.compareShapeletFunctions(sa, sb, rtolEllipse=rtolEllipse, rtolCoeff=rtolCoeff,
                                          atolEllipse=atolEllipse, atolCoeff=atolCoeff)

    def checkMoments(self, function, x, y, z):
        gx, gy = numpy.meshgrid(x, y)
        m = z.sum()
        dipole = lsst.afw.geom.Point2D((gx * z).sum() / m, (gy * z).sum() / m)
        gx -= dipole.getX()
        gy -= dipole.getY()
        quadrupole = lsst.afw.geom.ellipses.Quadrupole(
            (gx**2 * z).sum() / m,
            (gy**2 * z).sum() / m,
            (gx * gy * z).sum() / m
        )
        imageMoments = lsst.afw.geom.ellipses.Ellipse(quadrupole, dipole)
        shapeletMoments = function.evaluate().computeMoments()
        self.assertClose(imageMoments.getCenter().getX(), shapeletMoments.getCenter().getX(), rtol=1E-3)
        self.assertClose(imageMoments.getCenter().getY(), shapeletMoments.getCenter().getY(), rtol=1E-3)
        self.assertClose(imageMoments.getCore().getIxx(), shapeletMoments.getCore().getIxx(), rtol=1E-3)
        self.assertClose(imageMoments.getCore().getIyy(), shapeletMoments.getCore().getIyy(), rtol=1E-3)
        self.assertClose(imageMoments.getCore().getIxy(), shapeletMoments.getCore().getIxy(), rtol=1E-3)
        integral = numpy.trapz(numpy.trapz(z, gx, axis=1), y, axis=0)
        self.assertClose(integral, function.evaluate().integrate(), rtol=1E-3)

    def checkConvolution(self, f1, f2):
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-50, -50), lsst.afw.geom.Point2I(50, 50))
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
        if scipy is None:
            print "Skipping convolution test; scipy could not be imported."
            return
        # I'm using scipy.ndimage to convolve test images, because I can't figure
        # out how to make afw do it (afw can convolve images with kernels, but two similarly-sized
        # are apparently another matter; if I try to make a FixedKernel from one of the images,
        # I can't even make the operation commutative, let alone correct.
        scipy.ndimage.convolve(i1.getArray(), i2.getArray(), output=out.getArray(),
                               mode="constant", cval=0.0)
        self.assertClose(out.getArray(), ic1.getArray(), rtol=1E-4, atol=1E-5)
        self.assertClose(out.getArray(), ic2.getArray(), rtol=1E-4, atol=1E-5)
        return fc1, fc2
