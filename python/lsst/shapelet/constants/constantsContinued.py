from .._shapeletLib import BasisTypeEnum

from lsst.utils import continueClass

__all__ = []


@continueClass  # noqa: F811 (FIXME: remove for py 3.8+)
class BasisTypeEnum:  # noqa: F811
    # Workaround for broken pickling on Python 2
    # Without this fails with: TypeError: lsst.shapelet.constants.BasisTypeEnum.__new__(
    #     lsst.shapelet.constants.BasisTypeEnum) is not safe, use object.__new__()
    def __reduce__(self):
        return (BasisTypeEnum, (int(self), ))
