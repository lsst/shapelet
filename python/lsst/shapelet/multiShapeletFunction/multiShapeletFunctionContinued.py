from .multiShapeletFunction import MultiShapeletFunction

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa F811
class MultiShapeletFunction:
    def __reduce__(self):
        return (MultiShapeletFunction, (list(self.getComponents()),))
