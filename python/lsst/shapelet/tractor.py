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
import re
import warnings
try:
    import cPickle as pickle
except ImportError:
    import pickle

from .shapeletLib import *


def registerRadialProfiles():
    """Register the pickled profiles in the data directory with the RadialProfile singleton registry.

    This should only be called at import time by this module; it's only a function to avoid polluting
    the module namespace with all the local variables used here.
    """
    dataDir = os.path.join(os.environ["SHAPELET_DIR"], "data")
    regex = re.compile(r"([a-z]+\d?)_K(\d+)_MR(\d+)\.pickle")
    for filename in os.listdir(dataDir):
        match = regex.match(filename)
        if not match:
            continue
        name = match.group(1)
        nComponents = int(match.group(2))
        maxRadius = int(match.group(3))
        try:
            profile = RadialProfile.get(name)
        except lsst.pex.exceptions.Exception:
            warnings.warn("No C++ profile for multi-Gaussian pickle file '%s'" % filename)
            continue
        with open(os.path.join(dataDir, filename), 'r') as stream:
            array = pickle.load(stream)
        amplitudes = array[:nComponents]
        amplitudes /= amplitudes.sum()
        variances = array[nComponents:]
        if amplitudes.shape != (nComponents,) or variances.shape != (nComponents,):
            warnings.warn("Unknown format for multi-Gaussian pickle file '%s'" % filename)
            continue
        basis = MultiShapeletBasis(1)
        for amplitude, variance in zip(amplitudes, variances):
            radius = variance**0.5
            matrix = numpy.array([[amplitude / ShapeletFunction.FLUX_FACTOR]], dtype=float)
            basis.addComponent(radius, 0, matrix)
        profile.registerBasis(basis, nComponents, maxRadius)
# We register all the profiles at module import time, to allow C++ code to access all available profiles
# without having to later call Python code to unpickle them.
registerRadialProfiles()


def evaluateRadial(basis, r, sbNormalize=False, doComponents=False):
    """Plot a single-element MultiShapeletBasis as a radial profile.
    """
    ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes())
    coefficients = numpy.ones(1, dtype=float)
    msf = basis.makeFunction(ellipse, coefficients)
    ev = msf.evaluate()
    n = 1
    if doComponents:
        n += len(msf.getComponents())
    z = numpy.zeros((n,) + r.shape, dtype=float)
    for j, x in enumerate(r):
        z[0, j] = ev(x, 0.0)
    if doComponents:
        for i, sf in enumerate(msf.getComponents()):
            evc = sf.evaluate()
            for j, x in enumerate(r):
                z[i+1, j] = evc(x, 0.0)
    if sbNormalize:
        z /= ev(1.0, 0.0)
    return z


def integrateNormalizedFluxes(maxRadius=20.0, nSteps=5000):
    """!
    After normalizing by surface brightness at r=1 r_e, integrate the profiles to compare
    relative fluxes between the true profiles and their approximations.

    @param[in] maxRadius   Maximum radius to integrate the profile, in units of r_e.
    @param[in] nSteps      Number of concrete points at which to evaluate the profile to
                           do the integration (we just use the trapezoidal rule).
    """
    radii = numpy.linspace(0.0, maxRadius, nSteps)
    profiles = {name: RadialProfile.get(name) for name in ("exp", "lux", "dev", "luv",
                                                           "ser2", "ser3", "ser5")}
    evaluated = {}
    for name, profile in profiles.iteritems():
        evaluated[name] = profile.evaluate(radii)
        basis = profile.getBasis(8)
        evaluated["g" + name] = evaluateRadial(basis, radii, sbNormalize=True, doComponents=False)[0, :]
    fluxes = {name: numpy.trapz(z*radii, radii) for name, z in evaluated.iteritems()}
    return fluxes


def plotSuite(doComponents=False):
    """Plot all the profiles defined in this module together: true exp and dev, the SDSS softended/truncated
    lux and luv, and the multi-Gaussian approximations to all of these.

    To plot the individual Gaussians that form the multi-Gaussian approximations, pass doComponents=True.

    Returns a tuple of (figure, axes), where 'figure' is the matplotlib figure that contains the plot,
    and axes is a 2x4 NumPy array of matplotlib axes objects
    """
    from matplotlib import pyplot
    fig = pyplot.figure(figsize=(9, 4.7))
    axes = numpy.zeros((2, 4), dtype=object)
    r1 = numpy.logspace(-3, 0, 1000, base=10)
    r2 = numpy.linspace(1, 10, 1000)
    r = [r1, r2]
    for i in range(2):
        for j in range(4):
            axes[i, j] = fig.add_subplot(2, 4, i*4+j+1)
    profiles = {name: RadialProfile.get(name) for name in ("exp", "lux", "dev", "luv")}
    basis = {name: profiles[name].getBasis(8) for name in profiles}
    z = numpy.zeros((2, 4), dtype=object)
    colors = ("k", "g", "b", "r")
    fig.subplots_adjust(wspace=0.025, hspace=0.025, bottom=0.15, left=0.1, right=0.98, top=0.92)
    centers = [None, None]
    for i in range(2):   # 0=profile, 1=relative error
        for j in range(0, 4, 2):  # grid columns: 0=exp-like, 2=dev-like
            bbox0 = axes[i, j].get_position()
            bbox1 = axes[i, j+1].get_position()
            bbox1.x0 = bbox0.x1 - 0.06
            bbox0.x1 = bbox1.x0
            centers[j/2] = 0.5*(bbox0.x0 + bbox1.x1)
            axes[i, j].set_position(bbox0)
            axes[i, j+1].set_position(bbox1)
    for j in range(0, 2):
        z[0, j] = [evaluateRadial(basis[k], r[j], sbNormalize=True, doComponents=doComponents)
                   for k in ("exp", "lux")]
        z[0, j][0:0] = [profiles[k].evaluate(r[j])[numpy.newaxis, :] for k in ("exp", "lux")]
        z[0, j+2] = [evaluateRadial(basis[k], r[j], sbNormalize=True, doComponents=doComponents)
                     for k in ("dev", "luv")]
        z[0, j+2][0:0] = [profiles[k].evaluate(r[j])[numpy.newaxis, :] for k in ("dev", "luv")]
    methodNames = [["loglog", "semilogy"], ["semilogx", "plot"]]
    for j in range(0, 4):  # grid columns
        y[1, j] = [(y[0, j][0][0, :] - y[0, j][i][0, :])/y[0, j][0][0, :] for i in range(0, 4)]
        handles = []
        method0 = getattr(axes[0, j], methodNames[0][j%2])
        method1 = getattr(axes[1, j], methodNames[1][j%2])
        for k in range(4):
            y0 = y[0, j][k]
            handles.append(method0(r[j%2], y0[0, :], color=colors[k])[0])
            if doComponents:
                for l in range(1, y0.shape[0]):
                    method0(r[j%2], y0[l, :], color=colors[k], alpha=0.25)
            method1(r[j%2], y[1, j][k], color=colors[k])
        axes[0, j].set_xticklabels([])
        axes[0, j].set_ylim(1E-6, 1E3)
        axes[1, j].set_ylim(-0.2, 1.0)
    for i, label in enumerate(("profile", "relative error")):
        axes[i, 0].set_ylabel(label)
        for t in axes[i, 0].get_yticklabels():
            t.set_fontsize(11)
    for j in range(1, 4):
        axes[0, j].set_yticklabels([])
        axes[1, j].set_yticklabels([])
    xticks = [['$\\mathdefault{10^{%d}}$' % i for i in range(-3, 1)],
              [str(i) for i in range(1, 11)]]
    xticks[0][-1] = ""
    xticks[1][-1] = ""
    for j in range(0, 4):
        axes[1, j].set_xticklabels(xticks[j%2])
        for t in axes[1, j].get_xticklabels():
            t.set_fontsize(11)
    fig.legend(handles, ["exp/dev", "lux/luv", "approx exp/dev", "approx lux/luv"],
               loc='lower center', ncol=4)
    fig.text(centers[0], 0.95, "exponential", ha='center', weight='bold')
    fig.text(centers[1], 0.95, "de Vaucouleur", ha='center', weight='bold')
    return fig, axes
