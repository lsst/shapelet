from __future__ import absolute_import
from builtins import range
from builtins import object
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

from .shapeletLib import *
from . import tractor
from lsst.afw.geom.ellipses import Quadrupole as EllipseCore


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
