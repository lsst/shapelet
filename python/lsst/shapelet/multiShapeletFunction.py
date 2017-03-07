from __future__ import absolute_import

from ._multiShapeletFunction import MultiShapeletFunction

__all__ = []

def __reduce__(self):
    return (MultiShapeletFunction, (list(self.getComponents()),))
MultiShapeletFunction.__reduce__ = __reduce__
