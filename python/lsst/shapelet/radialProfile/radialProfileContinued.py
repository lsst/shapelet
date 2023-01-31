import numpy as np

from .._shapeletLib import RadialProfile

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa: F811 (FIXME: remove for py 3.8+)
class RadialProfile:  # noqa: F811
    def evaluate(self, r):
        if isinstance(r, np.ndarray):
            return self._evaluate(r.ravel()).reshape(r.shape)
        else:
            return self._evaluate(r)
