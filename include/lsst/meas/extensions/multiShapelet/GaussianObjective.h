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

protected:

    struct Component {
        double amplitude;
        PTR(afw::geom::ellipses::Ellipse) ellipse;
        Eigen::Matrix<double,5,Eigen::Dynamic> jacobian;

        explicit Component(int parameterSize) : amplitude(0.0), ellipse(), jacobian(5, parameterSize) {
            jacobian.setZero();
        }
    };

    typedef std::vector<Component> ComponentList;

    GaussianObjective(
        int nComponents, int parameterSize,
        afw::detection::Footprint const & region,
        ndarray::Array<double const,1,1> const & data,
        ndarray::Array<double const,1,1> const & weights = ndarray::Array<double,1,1>()
    );

    GaussianObjective(
        int nComponents, int parameterSize,
        afw::geom::Box2I const & bbox,
        ndarray::Array<double const,1,1> const & data,
        ndarray::Array<double const,1,1> const & weights = ndarray::Array<double,1,1>()
    );

    virtual void readParameters(
        ndarray::Array<double const,1,1> const & parameters, 
        ComponentList::iterator begin, ComponentList::iterator const end
    ) = 0;

private:

    typedef std::vector<GaussianModelBuilder> BuilderList;

    ComponentList _components;
    BuilderList _builders;
    ndarray::Array<double const,1,1> _data;
    ndarray::Array<double const,1,1> _weights;
};

}}}} // namespace lsst::meas::extensions::multisShapelet

#endif // !MULTISHAPELET_GaussianObjective_h_INCLUDED
