from __future__ import absolute_import

from ._constants import BasisTypeEnum

__all__ = []

# Workaround for broken pickling on Python 2
# Without this fails with: TypeError: lsst.shapelet._constants.BasisTypeEnum.__new__(
#     lsst.shapelet._constants.BasisTypeEnum) is not safe, use object.__new__()
def __reduce__(self):
    return (BasisTypeEnum, (int(self), ))
BasisTypeEnum.__reduce__ = __reduce__

