from __future__ import absolute_import, division, print_function

import numpy as np

from .radialProfile import RadialProfile

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa F811
class RadialProfile:
    def evaluate(self, r):
        if isinstance(r, np.ndarray):
            return self._evaluate(r.ravel()).reshape(r.shape)
        else:
            return self._evaluate(r)
