// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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
#ifndef MULTISHAPELET_GaussianObjective_h_INCLUDED
#define MULTISHAPELET_GaussianObjective_h_INCLUDED

#include "lsst/meas/extensions/multiShapelet/HybridOptimizer.h"
#include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class GaussianObjective : public Objective {
public:

    virtual void computeFunction(
        ndarray::Array<double const,1,1> const & parameters, 
        ndarray::Array<double,1,1> const & function
    );

    virtual void computeDerivative(
        ndarray::Array<double const,1,1> const & parameters, 
        ndarray::Array<double const,1,1> const & function,
        ndarray::Array<double,2,-2> const & derivative
    );

    double getAmplitude() const { return _amplitude; }
    
#ifndef SWIG // Don't need to create this class from Python, just use it.

    struct Component {
        double amplitude;
        double radius;

        explicit Component(double amplitude_, double radius_) : amplitude(amplitude_), radius(radius_) {}
    };

    typedef std::vector<Component> ComponentList;

    GaussianObjective(
        ComponentList const & components, afw::geom::Point2D const & center,
        afw::detection::Footprint const & region,
        ndarray::Array<double const,1,1> const & data,
        ndarray::Array<double const,1,1> const & weights = ndarray::Array<double,1,1>()
    );

    GaussianObjective(
        ComponentList const & components, afw::geom::Point2D const & center,
        afw::geom::Box2I const & bbox,
        ndarray::Array<double const,1,1> const & data,
        ndarray::Array<double const,1,1> const & weights = ndarray::Array<double,1,1>()
    );

#endif

private:

    typedef std::vector<GaussianModelBuilder> BuilderList;

    void _initialize();

    double _amplitude;
    double _modelSquaredNorm;
    afw::geom::ellipses::Ellipse _ellipse;
    ComponentList _components;
    BuilderList _builders;
    ndarray::Array<double,1,1> _model;
    ndarray::Array<double const,1,1> _data;
    ndarray::Array<double const,1,1> _weights;
};

}}}} // namespace lsst::meas::extensions::multisShapelet

#endif // !MULTISHAPELET_GaussianObjective_h_INCLUDED
