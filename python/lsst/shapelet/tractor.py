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

class sdss(object):
    """Namespace-only class for the softened/truncated exponential and de Vaucouleur
    profiles used in SDSS Photo model fits)."""

    DEFAC = -7.66925
    DEVOUT = 8.0
    DEVCUT = 7.0

    EXPFAC = -1.67835
    EXPOUT = 4.0
    EXPCUT = 3.0

    @classmethod
    def dev(cls, r):
        """Truncated de Vaucouleur - copied/translated from SDSS Photo package

        Expected input is a NumPy array of radius values.
        """
        p = numpy.exp(cls.DEFAC * ((r**2 + 0.0004)**0.125 - 1.0))
        big = r > cls.DEVCUT
        scr = (r[big] - cls.DEVCUT) / (cls.DEVOUT - cls.DEVCUT)
        scr = 1.0 - scr**2
        p[big] *= scr*scr
        p[r > cls.DEVOUT] = 0.0
        return p

    @classmethod
    def exp(cls, r):
        """Truncated exponential - copied/translated from SDSS from SDSS Photo package

        Expected input is a NumPy array of radius values.
        """
        p = numpy.exp(cls.EXPFAC * (r - 1.0))
        big = r > cls.EXPCUT
        scr = (r[big] - cls.EXPCUT) / (cls.EXPOUT - cls.EXPCUT);
        scr = 1.0 - scr**2
        p[big] *= scr * scr
        p[r > cls.EXPOUT] = 0.0
        return p

class exact(object):
    """Namespace-only class for exact exponential and de Vaucouleur profiles."""

    EXP_KAPPA = 1.67834699
    DEV_KAPPA = 7.66924944

    @classmethod
    def dev(cls, r):
        """de Vaucouleur (Sersic n=4) profile, normalized to unit surface brightness at r=1

        Expected input is a NumPy array of radius values.
        """
        return numpy.exp(-cls.DEV_KAPPA*(r**0.25 - 1.0))

    @classmethod
    def exp(cls, r):
        """exponential (Sersic n=1) profile, normalizedd to unit surface brightness at r=1

        Expected input is a NumPy array of radius values.
        """
        return numpy.exp(-cls.EXP_KAPPA*(r - 1.0))

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

def evaluateRadial(basis, r, sbNormalize=False, doComponents=False):
    """Return a radial profile for a single-element MultiShapeletBasis.
    """
    ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes())
    coefficients = numpy.ones(1, dtype=float)
    msf = basis.makeFunction(ellipse, coefficients)
    ev = msf.evaluate()
    n = 1
    if doComponents:
        n += len(msf.getElements())
    z = numpy.zeros((n,) + r.shape, dtype=float)
    for j, x in enumerate(r):
        z[0,j] = ev(x, 0.0)
    if doComponents:
        for i, sf in enumerate(msf.getElements()):
            evc = sf.evaluate()
            for j, x in enumerate(r):
                z[i+1,j] = evc(x, 0.0)
    if sbNormalize:
        z /= ev(1.0, 0.0)
    return z

def integrateNormalizedFluxes(maxRadius=20.0, nSteps=5000):
    radii = numpy.linspace(0.0, maxRadius, nSteps)
    basis = {name: loadBasis(name, nComponents=8) for name in ("exp", "lux", "dev", "luv")}
    profiles = {
        "exp": exact.exp(radii), "lux": sdss.exp(radii),
        "gexp": evaluateRadial(basis["exp"], radii, sbNormalize=True, doComponents=False)[0,:],
        "glux": evaluateRadial(basis["lux"], radii, sbNormalize=True, doComponents=False)[0,:],
        "dev": exact.dev(radii), "luv": sdss.dev(radii),
        "gdev": evaluateRadial(basis["dev"], radii, sbNormalize=True, doComponents=False)[0,:],
        "gluv": evaluateRadial(basis["luv"], radii, sbNormalize=True, doComponents=False)[0,:],
        }
    fluxes = {
        name: numpy.trapz(profile*radii, radii) for name, profile in profiles.iteritems()
        }
    return fluxes

def plotSuite(doComponents=False):
    from matplotlib import pyplot
    fig = pyplot.figure(figsize=(9,4.7))
    axes = numpy.zeros((2,4), dtype=object)
    r1 = numpy.logspace(-3, 0, 1000, base=10)
    r2 = numpy.linspace(1, 10, 1000)
    r = [r1, r2]
    for i in range(2):
        for j in range(4):
            axes[i,j] = fig.add_subplot(2, 4, i*4+j+1)
    basis = {name: loadBasis(name, nComponents=8) for name in ("exp", "lux", "dev", "luv")}
    z = numpy.zeros((2,4), dtype=object)
    colors = ("k", "g", "b", "r")
    fig.subplots_adjust(wspace=0.025, hspace=0.025, bottom=0.15, left=0.1, right=0.98, top=0.92)
    centers = [None, None]
    for i in range(2):
        for j in range(0,4,2):
            bbox0 = axes[i,j].get_position()
            bbox1 = axes[i,j+1].get_position()
            bbox1.x0 = bbox0.x1 - 0.06
            bbox0.x1 = bbox1.x0
            centers[j/2] = 0.5*(bbox0.x0 + bbox1.x1)
            axes[i,j].set_position(bbox0)
            axes[i,j+1].set_position(bbox1)
    for j in range(0,2):
        z[0,j] = [evaluateRadial(basis[k], r[j], sbNormalize=True, doComponents=doComponents)
                  for k in ("exp", "lux")]
        z[0,j][0:0] = [exact.exp(r[j])[numpy.newaxis,:], sdss.exp(r[j])[numpy.newaxis,:]]
        z[0,j+2] = [evaluateRadial(basis[k], r[j], sbNormalize=True, doComponents=doComponents)
                    for k in ("dev", "luv")]
        z[0,j+2][0:0] = [exact.dev(r[j])[numpy.newaxis,:], sdss.dev(r[j])[numpy.newaxis,:]]
    methodNames = [["loglog", "semilogy"], ["semilogx", "plot"]]
    for j in range(0,4):
        z[1,j] = [(z[0,j][0][0,:] - z[0,j][i][0,:])/z[0,j][0][0,:] for i in range(0,4)]
        handles = []
        method0 = getattr(axes[0,j], methodNames[0][j%2])
        method1 = getattr(axes[1,j], methodNames[1][j%2])
        for k in range(4):
            z0 = z[0,j][k]
            handles.append(method0(r[j%2], z0[0,:], color=colors[k])[0])
            if doComponents:
                for l in range(1, z0.shape[0]):
                    method0(r[j%2], z0[l,:], color=colors[k], alpha=0.25)
            method1(r[j%2], z[1,j][k], color=colors[k])
        axes[0,j].set_xticklabels([])
        axes[0,j].set_ylim(1E-6, 1E3)
        axes[1,j].set_ylim(-0.2, 1.0)
    for i, label in enumerate(("profile", "relative error")):
        axes[i,0].set_ylabel(label)
        for t in axes[i,0].get_yticklabels():
            t.set_fontsize(11)
    for j in range(1,4):
        axes[0,j].set_yticklabels([])
        axes[1,j].set_yticklabels([])
    xticks = [['$\\mathdefault{10^{%d}}$' % i for i in range(-3,1)],
               [str(i) for i in range(1,11)]]
    xticks[0][-1] = ""
    xticks[1][-1] = ""
    for j in range(0,4):
        axes[1,j].set_xticklabels(xticks[j%2])
        for t in axes[1,j].get_xticklabels():
            t.set_fontsize(11)
    fig.legend(handles, ["exp/dev", "lux/luv", "gexp/gdev", "glux/gluv"],
               loc='lower center', ncol=4)
    fig.text(centers[0], 0.95, "exponential", ha='center', weight='bold')
    fig.text(centers[1], 0.95, "de Vaucouleur", ha='center', weight='bold')
    return fig, axes
