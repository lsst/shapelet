from .._shapeletLib import MultiShapeletFunction

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa: F811 (FIXME: remove for py 3.8+)
class MultiShapeletFunction:  # noqa: F811
    def __reduce__(self):
        return (MultiShapeletFunction, (list(self.getComponents()),))
