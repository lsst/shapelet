from __future__ import absolute_import

import numpy as np

from ._radialProfile import RadialProfile

__all__ = []

def evaluate(self, r):
    if isinstance(r, np.ndarray):
        return self._evaluate(r.ravel()).reshape(r.shape)
    else:
        return self._evaluate(r)
RadialProfile.evaluate = evaluate

