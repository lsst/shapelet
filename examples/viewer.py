import lsst.meas.extensions.multiShapelet as ms
import lsst.afw.image
import numpy
from matplotlib import pyplot

class FitPsfViewer(object):

    class Iteration(object):

        def __init__(self, opt, viewer):
            self.state = opt.getState()
            if self.state & ms.HybridOptimizer.STEP_ACCEPTED:
                self.parameters = opt.getParameters().copy()
                self.chisq = opt.getChiSq()
            else:
                self.parameters = opt.getTrialParameters().copy()
                self.chisq = opt.getTrialChiSq()
            self.method = opt.getMethod()
            self.model = ms.FitPsfModel(viewer.ctrl, self.parameters)
            self.image = lsst.afw.image.ImageD(viewer.image.getBBox(lsst.afw.image.PARENT))
            self.model.evaluate(self.image, viewer.center)
            self.residuals = lsst.afw.image.ImageD(viewer.image, True)
            self.residuals -= self.image

        def __str__(self):
            return "state=%0x, chisq=%f, method=%s, parameters=%s" % (
                self.state, self.chisq, "LM" if self.method==ms.HybridOptimizer.LM else "BFGS",
                self.parameters
                )

        __repr__ = __str__

    def __init__(self, source, psf, ctrl=None):
        if ctrl is None:
            ctrl = ms.FitPsfControl()
        self.ctrl = ctrl
        self.saved = ms.FitPsfModel(ctrl, source)
        self.center = source.getCentroid()
        self.image = psf.computeImage(self.center)
        opt = ms.FitPsfAlgorithm.makeOptimizer(ctrl, self.image, self.center)
        maxIter = opt.getControl().maxIter
        self.iterations = [self.Iteration(opt, self)]
        for self.iterCount in range(maxIter):
            opt.step()
            self.iterations.append(self.Iteration(opt, self))
            if opt.getState() & ms.HybridOptimizer.FINISHED:
                break

    @staticmethod
    def _plotImage(image):
        bbox = image.getBBox(lsst.afw.image.PARENT)
        array = image.getArray()
        pyplot.imshow(array, interpolation='nearest', origin='lower',
                      extent=(bbox.getMinX()-0.5, bbox.getMaxX()+0.5, bbox.getMinY()-0.5, bbox.getMaxY()+0.5)
                      )
        ticks = [array.min(), 0.5 * (array.min() + array.max()), array.max()]
        pyplot.colorbar(orientation="horizontal", format="%.2g", ticks=ticks)        

    def plot(self, n):
        pyplot.clf()
        pyplot.subplot(1,3,1)
        self._plotImage(self.image)
        pyplot.subplot(1,3,2)
        self._plotImage(self.iterations[n].image)
        pyplot.subplot(1,3,3)
        self._plotImage(self.iterations[n].residuals)

