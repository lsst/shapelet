from __future__ import absolute_import

from ._shapeletFunction import ShapeletFunction

__all__ = []

def __reduce__(self):
    return (ShapeletFunction, (self.getOrder(), self.getBasisType(),
        self.getEllipse(), self.getCoefficients()))
ShapeletFunction.__reduce__ = __reduce__
