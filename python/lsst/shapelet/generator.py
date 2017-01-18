from __future__ import absolute_import
from builtins import range, object
from ._constants import HERMITE, LAGUERRE, computeSize

class IndexGenerator(object):
    """
    Base class for shapelet index generators.
    """

    __slots__ = "order", "size"

    @staticmethod
    def make(self, order, basisType):
        if basisType == HERMITE:
            return HermiteIndexGenerator(order)
        elif basisType == LAGUERRE:
            return LaguerreIndexGenerator(order)

    def __init__(self, order):
        self.order = order
        self.size = computeSize(self.order)

    def __len__(self):
        return self.size


class HermiteIndexGenerator(IndexGenerator):
    """
    Iterable that generates tuples of (i, nx, ny) in which:
     - 'i' is the overall coefficient index for a 2-d shapelet expansion (just counts from zero)
     - 'nx' is the order of the x expansion
     - 'ny' is the order of the y expansion
     """

    def __iter__(self):
        i = 0
        for n in range(0, self.order+1):
            for nx in range(0, n+1):
                yield (i, nx, n - nx)
                i += 1


class LaguerreIndexGenerator(IndexGenerator):
    """
    Iterable that generates tuples of (i, p, q, re) in which:
     - 'i' is the overall coefficient index for a 2-d shapelet expansion (just counts from zero)
     - 'p' and 'q' are the indices of the polar shapelet expansion (see BasisTypeEnum).
     - 're' is True if this the real part of the coefficient
     """

    def __iter__(self):
        i = 0
        for n in range(0, self.order+1):
            p = n
            q = 0
            while p > q:
                yield (i, p, q, True)
                i += 1
                yield (i, p, q, False)
                i += 1
                p -= 1
                q += 1
            if p == q:
                yield (i, p, q, True)
                i += 1
