// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "boost/make_shared.hpp"
#include "boost/math/special_functions/gamma.hpp"

#include "lsst/shapelet/RadialProfile.h"

namespace lsst { namespace shapelet {

namespace {

typedef std::map<std::string,RadialProfile*> RadialProfileRegistry;

RadialProfileRegistry & getRadialProfileRegistry() {
    static RadialProfileRegistry registry;
    return registry;
}


class SersicRadialProfile : public RadialProfile {
public:

    virtual double evaluate(double r) const {
        return std::exp(-_kappa*(std::pow(r, 1.0 / _n) - 1.0));
    }

    explicit SersicRadialProfile(std::string const & name, double n, int defaultMaxRadius=8) :
        RadialProfile(name, defaultMaxRadius), _n(n), _kappa(boost::math::gamma_q_inv(2*n, 0.5))
    {
        _momentsRadiusFactor = std::pow(_kappa, -n) * std::sqrt(boost::math::tgamma_ratio(4*n, 2*n)/2);
    }

private:
    double _n;
    double _kappa;
};

SersicRadialProfile expProfile("exp", 1.0);
SersicRadialProfile ser2Profile("ser2", 2.0);
SersicRadialProfile ser3Profile("ser3", 3.0);
SersicRadialProfile devProfile("dev", 4.0);
SersicRadialProfile ser5Profile("ser5", 5.0);

class GaussianRadialProfile : public SersicRadialProfile {
public:

    explicit GaussianRadialProfile() :
        SersicRadialProfile("gaussian", 0.5, 0)
    {
        PTR(MultiShapeletBasis) basis = boost::make_shared<MultiShapeletBasis>(1);
        ndarray::Array<double,2,2> matrix = ndarray::allocate(1,1);
        matrix.deep() = 0.5 / std::sqrt(afw::geom::PI);
        basis->addComponent(_momentsRadiusFactor, 0, matrix);
        registerBasis(basis, 1, 0);
    }

};

GaussianRadialProfile gaussianProfilex;

// Truncated exp profile used in SDSS pipeline
class LuxRadialProfile : public RadialProfile {
public:

    virtual double evaluate(double r) const {
        // formula and constants taken from SDSS Photo code
        static double const EXPFAC = -1.67835;
        static double const EXPOUT = 4.0;
        static double const EXPCUT = 3.0;
        if (r > EXPOUT) return 0.0;
        double p = std::exp(EXPFAC * (r - 1.0));
        if (r > EXPCUT) {
            double s = (r - EXPCUT) / (EXPOUT - EXPCUT);
            s = 1.0 - s*s;
            p *= s*s;
        }
        return p;
    }

    LuxRadialProfile() : RadialProfile("lux", 4) {
        _momentsRadiusFactor = 0.95473; // computed in Mathematica
    }
};

LuxRadialProfile luxProfile;

// Softened, truncated de Vaucouleur profile used in SDSS pipeline
class LuvRadialProfile : public RadialProfile {
public:

    virtual double evaluate(double r) const {
        // formula and constants taken from SDSS Photo code
        static double const DEFAC = -7.66925;
        static double const DEVOUT = 8.0;
        static double const DEVCUT = 7.0;
        if (r > DEVOUT) return 0.0;
        double p = std::exp(DEFAC*(std::pow(r*r + 0.0004, 0.125) - 1.0));
        if (r > DEVCUT) {
            double s = (r - DEVCUT) / (DEVOUT - DEVCUT);
            s = 1.0 - s*s;
            p *= s*s;
        }
        return p;
    }

    LuvRadialProfile() : RadialProfile("luv", 8) {
        _momentsRadiusFactor = 1.57016; // computed in Mathematica
    }
};

LuvRadialProfile luvProfile;

} // anonymous

RadialProfile::RadialProfile(std::string const & name, int defaultMaxRadius) :
    _momentsRadiusFactor(std::numeric_limits<double>::quiet_NaN()),
    _name(name),
    _defaultMaxRadius(defaultMaxRadius)
{
    getRadialProfileRegistry()[_name] = this;
}

RadialProfile & RadialProfile::get(std::string const & name) {
    RadialProfileRegistry::const_iterator i = getRadialProfileRegistry().find(name);
    if (i == getRadialProfileRegistry().end()) {
        throw LSST_EXCEPT(
            pex::exceptions::NotFoundException,
            (boost::format("RadialProfile with name '%s' not found") % name).str()
        );
    }
    return *i->second;
}

ndarray::Array<double,1,1> RadialProfile::evaluate(ndarray::Array<double const,1,1> const & r) const {
    ndarray::Array<double,1,1> p = ndarray::allocate(r.getSize<0>());
    for (int i = 0, n = r.getSize<0>(); i < n; ++i) {
        p[i] = evaluate(r[i]);
    }
    return p;
}

PTR(MultiShapeletBasis) RadialProfile::getBasis(int nComponents, int maxRadius) const {
    if (maxRadius == 0) maxRadius = _defaultMaxRadius;
    BasisRegistry::const_iterator i = _basisRegistry.find(std::make_pair(nComponents, maxRadius));
    if (i == _basisRegistry.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::NotFoundException,
            (boost::format("Basis for profile '%s' with nComponents=%d and maxRadius=%d not present")
             % _name % nComponents % maxRadius).str()
        );
    }
    // Ideally, we'd return a const ptr instead of copying, but because Swig casts away the
    // constness, that's not sufficiently safe.
    return boost::make_shared<MultiShapeletBasis>(*i->second);
}

void RadialProfile::registerBasis(PTR(MultiShapeletBasis) basis, int nComponents, int maxRadius) {
    _basisRegistry[std::make_pair(nComponents, maxRadius)] = basis;
}



}} // namespace lsst::shapelet
