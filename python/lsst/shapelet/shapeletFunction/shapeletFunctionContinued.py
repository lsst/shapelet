from .._shapeletLib import ShapeletFunction

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa: F811 (FIXME: remove for py 3.8+)
class ShapeletFunction:  # noqa: F811
    def __reduce__(self):
        return (ShapeletFunction, (self.getOrder(), self.getBasisType(),
                                   self.getEllipse(), self.getCoefficients()))
