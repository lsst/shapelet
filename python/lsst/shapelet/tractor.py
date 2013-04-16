#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
"""
Code to load multi-Gaussian approximations to profiles from "The Tractor"
into a lsst.shapelet.MultiShapeletBasis.

Please see the README file in the data directory of the lsst.shapelet
package for more information.
"""

import numpy
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

import lsst.shapelet

def loadParameters(profile, nComponents, maxRadius=None):
    """Load the parameters of a multi-Gaussian profile

    @param[in]   profile       Name of the ideal profile being approximated; usually one of
                               'exp', 'dev', 'ser2', 'ser3', 'lux', or 'luv'.
    @param[in]   nComponents   Number of Gaussians in the approximation.
    @param[in]   maxRadius     Maximum radius at which the multi-Gaussian approximated was optimized, in
                               units of the half-light radius.  Defaults to 8 for all profiles except 'lux',
                               which is only defined to 4.

    @return a tuple of (amplitudes, variances), each of which are NumPy arrays of length nComponents,
            containing the relative flux and variance of the components

    The returned profiles are defined such that the half-light radius of the profile and integrated flux
    are both unity.
    """
    if maxRadius is None:
        if profile == 'lux':
            maxRadius = 4
        else:
            maxRadius = 8
    name = "%s_K%02d_MR%02d.pickle" % (profile, nComponents, maxRadius)
    path = os.path.join(os.environ["SHAPELET_DIR"], "data", name)
    with open(path, 'r') as stream:
        array = pickle.load(stream)
    amplitudes = array[:nComponents]
    amplitudes /= amplitudes.sum()
    variances = array[nComponents:]
    if amplitudes.shape != (nComponents,) or variances.shape != (nComponents,):
        raise ValueError("Unknown format for pickled profile file %s" % name)
    return amplitudes, variances

def loadBasis(profile, nComponents, maxRadius=None):
    """Load the parameters of a multi-Gaussian profile into a single-element MultiShapeletBasis.

    @param[in]   profile       Name of the ideal profile being approximated; usually one of
                               'exp', 'dev', 'ser2', 'ser3', 'lux', or 'luv'.
    @param[in]   nComponents   Number of Gaussians in the approximation.
    @param[in]   maxRadius     Maximum radius at which the multi-Gaussian approximated was optimized, in
                               units of the half-light radius.  Defaults to 8 for all profiles except 'lux',
                               which is only defined to 4.

    @return a MultiShapeletBasis object

    The returned profiles are defined such that the half-light radius of the profile and integrated flux
    are both unity.
    """
    amplitudes, variances = loadParameters(profile, nComponents, maxRadius)
    basis = lsst.shapelet.MultiShapeletBasis(1)
    for amplitude, variance in zip(amplitudes, variances):
        radius = variance**0.5
        matrix = numpy.array([[amplitude / lsst.shapelet.ShapeletFunction.FLUX_FACTOR]], dtype=float)
        basis.addComponent(radius, 0, matrix)
    return basis

def plotRadial(basis):
    """Plot a single-element MultiShapeletBasis as a radial profile.
    """
    from matplotlib import pyplot
    ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes())
    coefficients = numpy.ones(1, dtype=float)
    msf = basis.makeFunction(ellipse, coefficients)
    ev = msf.evaluate()
    r = numpy.linspace(0, 8, 1000)
    z = numpy.zeros(r.shape, dtype=float)
    for n, x in enumerate(r):
        z[n] = ev(x, 0.0)
    pyplot.plot(r, z)
    return r, z
