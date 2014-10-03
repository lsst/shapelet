// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_SHAPELET_RadialProfile_h_INCLUDED
#define LSST_SHAPELET_RadialProfile_h_INCLUDED

#include <map>

#include "lsst/shapelet/MultiShapeletBasis.h"

namespace lsst { namespace shapelet {

/**
 *  Registry and utility class for multi-Gaussian approximations to radial profiles.
 *
 *  RadialProfile provides C++ access to the saved multi-Gaussian approximations stored in the pickle
 *  files in the data subdirectory of the shapelet package, by maintaining a singleton registry of these
 *  that is populated at import time by tractor.py.  It also provides access to the true (exact) profiles,
 *  (mostly useful for educational/diagnostic purposes) and information about the relationship between
 *  different radii for this profile (see getMomentsRadiusFactor()).
 *
 *  Each RadialProfile object represents a particular "true" profile, and can hold multiple multi-Gaussian
 *  (or multi-Shapelet) approximations with different degrees of fidelity, as controlled by the number of
 *  components in the approximation and the maximum radius at which the approximation was fit.  For more
 *  information about these approximations, see tractor.py and references therein.  All RadialProfiles
 *  are defined in units of the half-light radius of the true profile, even for profiles for which
 *  this is not the radius typically used.
 *
 *  Several predefined subclasses of RadialProfile are defined privately, and can be accessed by name
 *  through the static get() method:
 *   - "gaussian": a single Gaussian profile
 *   - "exp": Exponential profile (Sersic n=1)
 *   - "ser2": Sersic profile with n=2
 *   - "ser3": Sersic profile with n=3
 *   - "dev": de Vaucouleur profile (Sersic n=4)
 *   - "ser5": Sersic profile with n=5
 *   - "lux": SDSS truncated Exponential profile
 *   - "luv": SDSS truncated and softened de Vaucouleur profile
 */
class RadialProfile {
public:

    /**
     *  @brief Return a predefined radial profile given its name.
     *
     *  New RadialProfile classes are registered when their constructor is invoked, so all subclasses
     *  should only be instantiated once (at module scope).
     */
    static RadialProfile & get(std::string const & name);

    /// Return the name of the profile.
    std::string const getName() const { return _name; }

    /**
     *  @brief Evaluate the profile at the given radius
     *
     *  Radius is given in units of half-light radius, and the profile is normalized to 1 at r=1.
     */
    virtual double evaluate(double r) const = 0;

    /**
     *  @brief Evaluate the profile at the given radius (vectorized).
     *
     *  Radius is given in units of half-light radius, and the profile is normalized to 1 at r=1.
     */
    ndarray::Array<double,1,1> evaluate(ndarray::Array<double const,1,1> const & r) const;

    /// Return the 2nd-moment radius in units of the half-light radius.
    double getMomentsRadiusFactor() const { return _momentsRadiusFactor; }

    /**
     *  @brief Return a multi-Gaussian approximation to the radial profile
     *
     *  Unlike the evaluate() method, the returned basis is normalized that the integrated
     *  flux is equal to the (single) amplitude.
     *
     *  @param[in] nComponents   Number of Gaussians used to approximate the profile.
     *  @param[in] maxRadius     Maximum radius used in fitting the approximation.  Passing 0
     *                           selects the default for that profile.
     *
     *  Gaussian approximations are computed in advance and persisted, and are usually loaded
     *  and registered with the RadialProfile object in Python at module-import time.
     *
     *  Throws NotFoundError if the a basis with the given combination of nComponents and
     *  maxRadius has not been added to the RadialProfile.
     */
    PTR(MultiShapeletBasis) getBasis(int nComponents, int maxRadius=0) const;

    /**
     *  @brief Register a basis with the profile, making it available to callers of getBasis().
     */
    void registerBasis(PTR(MultiShapeletBasis) basis, int nComponents, int maxRadius);

protected:

    explicit RadialProfile(std::string const & name, int defaultMaxRadius);

    double _momentsRadiusFactor;

private:

    typedef std::map<std::pair<int,int>,PTR(MultiShapeletBasis)> BasisRegistry;

    std::string const _name;
    int _defaultMaxRadius;
    BasisRegistry _basisRegistry;
};

}} // namespace lsst::shapelet

#endif // !LSST_SHAPELET_RadialProfile_h_INCLUDED
