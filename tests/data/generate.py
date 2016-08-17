#!/usr/bin/env python
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
Script to generate test data.

Test data is checked into git rather than generated on the fly
because generating it requires the GalSim package
(https://github.com/GalSim-developers/GalSim), which is not part
of the LSST stack.  But it's still better to save the script as
well, even if it can't usually be run.
"""

import galsim
import lsst.shapelet.tractor

GALAXY_RADIUS = 8.0
PSF_SIGMA = 2.0
E1 = 0.3
E2 = -0.2
PROFILES = [
    ("exp", galsim.Exponential(half_light_radius=GALAXY_RADIUS), 9, 8),
    ("ser2", galsim.Sersic(n=2.0, half_light_radius=GALAXY_RADIUS), 9, 8),
    ("ser3", galsim.Sersic(n=3.0, half_light_radius=GALAXY_RADIUS), 9, 8),
    ("dev", galsim.DeVaucouleurs(half_light_radius=GALAXY_RADIUS), 9, 8),
]

MAKE_COMPONENT_IMAGES = False


def makeGaussian(flux, e1, e2, sigma):
    g = galsim.Gaussian(flux=flux, sigma=sigma)
    g.applyShear(e1=e1, e2=e2)
    return g


def main():
    for profile, exactGal, nComponents, maxRadius in PROFILES:
        exactGal.applyShear(e1=E1, e2=E2)
        psf = galsim.Gaussian(sigma=PSF_SIGMA)
        exactObj = galsim.Convolve([exactGal, psf])
        exactImage = exactObj.draw(dx=1.0)
        exactImage.write(profile + "-exact.fits")
        amplitudes, variances = lsst.shapelet.tractor.loadParameters(profile, nComponents, maxRadius)
        componentGals = []
        for n, (amplitude, variance) in enumerate(zip(amplitudes, variances)):
            componentGal = galsim.Gaussian(flux=amplitude, sigma=GALAXY_RADIUS*variance**0.5)
            componentGal.applyShear(e1=E1, e2=E2)
            componentGals.append(componentGal)
            if MAKE_COMPONENT_IMAGES:
                componentObj = galsim.Convolve([componentGal, psf])
                componentImage = galsim.ImageD(exactImage.bounds)
                componentObj.draw(componentImage, dx=1.0)
                componentImage.write("%s-approx-%0d.fits" % (profile, n))
        approxGal = galsim.Add(componentGals)
        approxObj = galsim.Convolve([approxGal, psf])
        approxImage = galsim.ImageD(exactImage.bounds)
        approxObj.draw(approxImage, dx=1.0)
        approxImage.write(profile + "-approx.fits")

    f1 = galsim.Add([
        makeGaussian(flux=0.600000, e1=0.09090909, e2=0.36363636, sigma=2.25810086),
        makeGaussian(flux=0.400000, e1=-0.11111111, e2=-0.11111111, sigma=2.98130750),
    ])
    f2 = galsim.Add([
        makeGaussian(flux=0.350000, e1=-0.26315789, e2=-0.21052632, sigma=2.99069756),
        makeGaussian(flux=0.650000, e1=-0.12500000, e2=0.12500000, sigma=2.80606626),
    ])
    f3 = galsim.Convolve([f1, f2])
    image = galsim.ImageD(galsim.BoundsI(-20, 20, -20, 20))
    f3.draw(image, dx=1.0)
    image.write("gaussians.fits")

if __name__ == "__main__":
    main()
