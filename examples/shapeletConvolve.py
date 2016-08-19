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


def plotShapeletFunction(axes, func, x, y):
    z = numpy.zeros((y.size, x.size), dtype=float)
    ev = func.evaluate()
    for i, py in enumerate(y):
        for j, px in enumerate(x):
            z[i, j] = ev(float(px), float(py))
    axes.imshow(z, interpolation='nearest', origin='lower')
    axes.yaxis.set_major_locator(ticker.NullLocator())
    axes.xaxis.set_major_locator(ticker.NullLocator())


def main():
    x = numpy.linspace(-5, 5, 101)
    y = numpy.linspace(-5, 5, 101)
    ellipse1 = ellipses.Ellipse(ellipses.Axes(1.0, 1.0, 0.3))
    ellipse2 = ellipses.Ellipse(ellipses.Axes(1.0, 1.0, numpy.pi/2 + 0.3))
    f1 = lsst.shapelet.ShapeletFunction(1, lsst.shapelet.HERMITE)
    f1.getCoefficients()[1] = 1.0
    f1.setEllipse(ellipse1)
    f2 = lsst.shapelet.ShapeletFunction(2, lsst.shapelet.HERMITE)
    f2.getCoefficients()[4] = 1.0
    f2.setEllipse(ellipse2)
    fC = f1.convolve(f2)
    plotShapeletFunction(pyplot.subplot(1, 3, 1), f1, x, y)
    plotShapeletFunction(pyplot.subplot(1, 3, 2), f2, x, y)
    plotShapeletFunction(pyplot.subplot(1, 3, 3), fC, x, y)
    pyplot.show()

if __name__ == "__main__":
    numpy.set_printoptions(suppress=True)
    main()
