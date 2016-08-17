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

import lsst.shapelet
from lsst.afw import geom
from lsst.afw.geom import ellipses
from matplotlib import pyplot
from matplotlib import ticker
import numpy


def makeBasisImages(basis, x, y):
    z = numpy.zeros((y.size, x.size, lsst.shapelet.computeSize(basis.getOrder())), dtype=float)
    for i, py in enumerate(y):
        for j, px in enumerate(x):
            basis.fillEvaluation(z[i, j, :], float(px), float(py))
    return z


def compareMoments(basis, x, y, z):
    e = basis.evaluate()
    monopole = z.sum()
    dipole = geom.Point2D((x * z).sum() / monopole, (y * z).sum() / monopole)
    dx = x - dipole.getX()
    dy = y - dipole.getY()
    quadrupole = ellipses.Quadrupole(
        (dx**2 * z).sum() / monopole,
        (dy**1 * z).sum() / monopole,
        (dx * dy * z).sum() / monopole
    )
    print ellipses.Ellipse(quadrupole, monopole)
    print e.computeMoments()


def plotBasisImages(basis, z):
    n = basis.getOrder()
    k = 0
    vmin = z.min()
    vmax = z.max()
    pyplot.figure()
    for i in range(n+1):
        for j in range(i+1):
            axes = pyplot.subplot(n+1, n+1, (n+1) * i + j + 1)
            axes.imshow(z[:, :, k], interpolation='nearest', origin='lower', vmin=vmin, vmax=vmax)
            axes.yaxis.set_major_locator(ticker.NullLocator())
            axes.xaxis.set_major_locator(ticker.NullLocator())
            if basis.getBasisType() == lsst.shapelet.HERMITE:
                pyplot.xlabel("x=%d, y=%d" % (j, i-j))
            else:
                pyplot.xlabel("p=%d, q=%d (%s)" % (i-j/2, j/2, "Im" if j % 2 else "Re"))
            k += 1


def processBasis(basis, x, y):
    z = makeBasisImages(basis, x, y)
    plotBasisImages(basis, z)


def main():
    x = numpy.linspace(-5, 5, 101)
    y = numpy.linspace(-5, 5, 101)
    ellipse = ellipses.Quadrupole(ellipses.Axes(1.2, 0.8, 0.3))
    hermiteBasis = lsst.shapelet.BasisEvaluator(4, lsst.shapelet.HERMITE)
    laguerreBasis = lsst.shapelet.BasisEvaluator(4, lsst.shapelet.LAGUERRE)
    processBasis(hermiteBasis, x, y)
    processBasis(laguerreBasis, x, y)
    pyplot.show()

if __name__ == "__main__":
    main()
