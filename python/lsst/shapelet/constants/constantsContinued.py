from __future__ import absolute_import, division, print_function

from .constants import BasisTypeEnum

from lsst.utils import continueClass

__all__ = []

@continueClass
class BasisTypeEnum:
    # Workaround for broken pickling on Python 2
    # Without this fails with: TypeError: lsst.shapelet.constants.BasisTypeEnum.__new__(
    #     lsst.shapelet.constants.BasisTypeEnum) is not safe, use object.__new__()
    def __reduce__(self):
        return (BasisTypeEnum, (int(self), ))

