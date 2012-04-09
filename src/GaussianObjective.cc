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

#include "lsst/meas/extensions/multiShapelet/GaussianObjective.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

void GaussianObjective::computeFunction(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double,1,1> const & function
) {
    readParameters(parameters, _components.begin(), _components.end());
    ndarray::EigenView<double,1,1> f(function);
    f.setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _builders[n].update(_components[n].ellipse);
        f -= _components[n].amplitude * _builders[n].getModel().asEigen();
    }
    f += _data.asEigen();
    if (!_weights.isEmpty()) {
        f.array() *= _weights.asEigen<Eigen::ArrayXpr>();
    }
}

void GaussianObjective::computeDerivative(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double const,1,1> const & function,
    ndarray::Array<double,2,-2> const & derivative
) {
    derivative.asEigen().setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _components[n].jacobian *= _components[n].amplitude;
        _builders[n].computeDerivative(derivative, _components[n].jacobian, true);
    }
}

GaussianObjective::GaussianObjective(
    int nComponents, int parameterSize,
    afw::detection::Footprint const & region,
    ndarray::Array<double const,1,1> const & data,
    ndarray::Array<double const,1,1> const & weights
) : Objective(region.getArea(), parameterSize),
    _components(nComponents), _builders(nComponents, GaussianModelBuilder(region)), 
    _data(data), _weights(weights)
{
    assert(getFunctionSize() == _data.getSize<0>());
    assert(_weights.isEmpty() || getFunctionSize() == _weights.getSize<0>());
}

GaussianObjective::GaussianObjective(
    int nComponents, int parameterSize,
    afw::geom::Box2I const & bbox,
    ndarray::Array<double const,1,1> const & data,
    ndarray::Array<double const,1,1> const & weights
) : Objective(bbox.getArea(), parameterSize),
    _components(nComponents), _builders(nComponents, GaussianModelBuilder(bbox)),
    _data(data), _weights(weights)
{
    assert(getFunctionSize() == _data.getSize<0>());
    assert(_weights.isEmpty() || getFunctionSize() == _weights.getSize<0>());
}

}}}} // namespace lsst::meas::extensions::multiShapelet
