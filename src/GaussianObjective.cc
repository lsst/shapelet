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
    _ellipse.getCore().readParameters(parameters.getData());
    ndarray::EigenView<double,1,1> model(_model);
    model.setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _ellipse.getCore().scale(_components[n].radius);
        model += _components[n].amplitude * _builders[n].computeModel(_ellipse).asEigen();
        _ellipse.getCore().scale(1.0 / _components[n].radius);
    }
    if (!_weights.isEmpty()) {
        model.array() *= _weights.asEigen<Eigen::ArrayXpr>();
    }
    _modelSquaredNorm = model.squaredNorm();
    _amplitude = model.dot(_data.asEigen()) / _modelSquaredNorm;
    function.asEigen() = _amplitude * model - _data.asEigen();
}

void GaussianObjective::computeDerivative(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double const,1,1> const & function,
    ndarray::Array<double,2,-2> const & derivative
) {
    derivative.asEigen().setZero();
    Eigen::Matrix<double,5,Eigen::Dynamic> jacobian(5, 3);
    jacobian.setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        afw::geom::LinearTransform scaling = afw::geom::LinearTransform::makeScaling(_components[n].radius);
        jacobian.block<3,3>(0, 0) = _ellipse.getCore().transform(scaling).d() * _components[n].amplitude;
        _ellipse.getCore().scale(_components[n].radius);
        _builders[n].computeDerivative(derivative, _ellipse, jacobian, true, false);
        _ellipse.getCore().scale(1.0 /_components[n].radius);
    }
    if (!_weights.isEmpty()) {
        derivative.asEigen<Eigen::ArrayXpr>() 
            *= (_weights.asEigen() * Eigen::RowVectorXd::Ones(parameters.getSize<0>())).array();
    }
    // Right now, 'derivative' is the partial derivative w.r.t. the objective parameters
    // with amplitude held fixed at 1.  However, the parameters also affect the amplitude, so we need
    // to compute the partial derivative of that.
    Eigen::VectorXd tmp = _data.asEigen() - 2.0 * _amplitude * _model.asEigen();
    Eigen::VectorXd dAmplitude = (derivative.asEigen().adjoint() * tmp) / _modelSquaredNorm;
    // Now we update 'derivative' so it becomes the complete derivative rather than the partial.
    derivative.asEigen() *= _amplitude;
    derivative.asEigen() += _model.asEigen() * dAmplitude.transpose();
}

GaussianObjective::GaussianObjective(
    ComponentList const & components, afw::geom::Point2D const & center,
    afw::detection::Footprint const & region,
    ndarray::Array<double const,1,1> const & data,
    ndarray::Array<double const,1,1> const & weights
) : Objective(region.getArea(), 3), _amplitude(1.0), _modelSquaredNorm(1.0),
    _ellipse(afw::geom::ellipses::Axes(), center),
    _components(components), _builders(components.size(), GaussianModelBuilder(region)), 
    _model(ndarray::allocate(region.getArea())), _data(data), _weights(weights)
{
    _initialize();
}

GaussianObjective::GaussianObjective(
    ComponentList const & components, afw::geom::Point2D const & center,
    afw::geom::Box2I const & bbox,
    ndarray::Array<double const,1,1> const & data,
    ndarray::Array<double const,1,1> const & weights
) : Objective(bbox.getArea(), 3), _amplitude(1.0), _modelSquaredNorm(1.0),
    _ellipse(afw::geom::ellipses::Axes(), center),
    _components(components), _builders(components.size(), GaussianModelBuilder(bbox)),
    _model(ndarray::allocate(bbox.getArea())), _data(data), _weights(weights)
{
    _initialize();
}

void GaussianObjective::_initialize() {
    assert(getFunctionSize() == _data.getSize<0>());
    assert(_weights.isEmpty() || getFunctionSize() == _weights.getSize<0>());
    if (!_weights.isEmpty()) {
        ndarray::EigenView<double,1,1,Eigen::ArrayXpr> wd(ndarray::copy(_data));
        wd *= _weights.asEigen<Eigen::ArrayXpr>();
        _data = wd.shallow();
    }
    _model.deep() = 0.0;
}

}}}} // namespace lsst::meas::extensions::multiShapelet
