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
            z[i,j] = ev(float(px), float(py))
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
