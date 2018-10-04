from .shapeletFunction import ShapeletFunction

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa
class ShapeletFunction:
    def __reduce__(self):
        return (ShapeletFunction, (self.getOrder(), self.getBasisType(),
                                   self.getEllipse(), self.getCoefficients()))
