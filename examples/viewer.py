import lsst.meas.extensions.multiShapelet as ms
import collections

class FitPsfViewer(object):

    Iteration = collections.namedtuple("Iteration", "model state chisq", verbose=False)

    def __init__(self, source, psf, ctrl=None):
        if ctrl is None:
            ctrl = ms.FitPsfControl()
        self.saved = ms.FitPsfModel(ctrl, source)
        self.center = source.getCentroid()
        self.image = psf.computeImage(self.center)
        opt = ms.FitPsfAlgorithm.makeOptimizer(ctrl, self.image, self.center)
        self.initial = ms.FitPsfModel(ctrl, opt.getParameters())
        maxIter = opt.getControl().maxIter
        self.iterations = []
        for self.iterCount in range(maxIter):
            opt.step()
            model = ms.FitPsfModel(ctrl, opt.getParameters())
            self.iterations.append(self.Iteration(model, opt.getState(), opt.getChiSq()))
            if opt.getState():
                break
        self.final = self.iterations[-1].model
