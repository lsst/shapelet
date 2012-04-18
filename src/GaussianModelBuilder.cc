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

#include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

GaussianModelBuilder::GaussianModelBuilder(afw::detection::Footprint const & region) :
    _xy(region.getArea(), 2), _xyt(region.getArea(), 2)
{
    int n = 0;
    for (
        afw::detection::Footprint::SpanList::const_iterator i = region.getSpans().begin();
        i != region.getSpans().end();
        ++i
    ) {
        for (int x = (**i).getX0(); x <= (**i).getX1(); ++x, ++n) {
            _xy(n, 0) = x;
            _xy(n, 1) = (**i).getY();
        }
    }
}

GaussianModelBuilder::GaussianModelBuilder(afw::geom::Box2I const & region) :
    _xy(region.getArea(), 2), _xyt(region.getArea(), 2)
{
    int n = 0;
    afw::geom::Point2I const llc = region.getMin();
    afw::geom::Point2I const urc = region.getMax();
    for (int y = llc.getY(); y <= urc.getY(); ++y) {
        for (int x = llc.getX(); x <= urc.getX(); ++x, ++n) {
            _xy(n, 0) = x;
            _xy(n, 1) = y;
        }
    }
}

GaussianModelBuilder::GaussianModelBuilder(GaussianModelBuilder const & other) :
    _ellipse(other._ellipse), _model(other._model), _xy(other._xy), _xyt(other._xyt)
{
    if (_ellipse) _ellipse = boost::make_shared<afw::geom::ellipses::Ellipse>(*_ellipse);
    if (!_model.isEmpty()) _model = ndarray::copy(_model);
}

GaussianModelBuilder & GaussianModelBuilder::operator=(GaussianModelBuilder const & other) {
    if (&other != this) {
        _ellipse = other._ellipse;
        _model = other._model;
        _xy = other._xy;
        _xyt = other._xyt;
        if (_ellipse) _ellipse = boost::make_shared<afw::geom::ellipses::Ellipse>(*_ellipse);
        if (!_model.isEmpty()) _model = ndarray::copy(_model);
    }
    return *this;
}

void GaussianModelBuilder::update(PTR(afw::geom::ellipses::Ellipse) const & ellipse) {
    _ellipse = ellipse;
    assert(_ellipse);
    if (_model.isEmpty()) {
        _model = ndarray::allocate(_xy.rows());
    }
    afw::geom::AffineTransform transform = _ellipse->getGridTransform();
    Eigen::Matrix2d const m = transform.getLinear().getMatrix();
    Eigen::Vector2d const t = transform.getTranslation().asEigen();
    _xyt.transpose() =  m * _xy.transpose();
    _xyt.transpose().colwise() += t;
    _model.asEigen<Eigen::ArrayXpr>() = std::exp(-0.5 * _xyt.rowwise().squaredNorm().array());
}

void GaussianModelBuilder::update(afw::geom::ellipses::Ellipse const & ellipse) {
    update(boost::make_shared<afw::geom::ellipses::Ellipse>(ellipse));
}

void GaussianModelBuilder::computeDerivative(ndarray::Array<double,2,-1> const & output) const {
    if (!_ellipse) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "update() must be called before computeDerivative"
        );
    }
    Eigen::Matrix<double,6,Eigen::Dynamic> gtJac(6, 5);
    gtJac.block<6,5>(0,0) = _ellipse->getGridTransform().d();
    _computeDerivative(output, gtJac, false);
}

void GaussianModelBuilder::computeDerivative(
    ndarray::Array<double,2,-1> const & output,
    Eigen::Matrix<double,5,Eigen::Dynamic> const & jacobian,
    bool add
) const {
    if (!_ellipse) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "update() must be called before computeDerivative"
        );
    }
    afw::geom::ellipses::Ellipse::GridTransform::DerivativeMatrix gtJac = _ellipse->getGridTransform().d();
    Eigen::Matrix<double,6,Eigen::Dynamic> finalJac = gtJac * jacobian;
    _computeDerivative(output, finalJac, add);
}

void GaussianModelBuilder::_computeDerivative(
    ndarray::Array<double,2,-1> const & output,
    Eigen::Matrix<double,6,Eigen::Dynamic> const & jacobian,
    bool add
) const {
    typedef afw::geom::AffineTransform AT;
    if (output.getSize<0>() != _xy.rows()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect number of rows for array: got %d, expected %d")
             % output.getSize<0>() % _xy.rows()).str()
        );
    }
    if (output.getSize<0>() != _xy.rows()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Mismatch between array (%d) and jacobian dimensions (%d)")
             % output.getSize<1>() % jacobian.cols()).str()
        );
    }
    Eigen::ArrayXd dfdx = -_xyt.col(0).array() * _model.asEigen<Eigen::ArrayXpr>();
    Eigen::ArrayXd dfdy = -_xyt.col(1).array() * _model.asEigen<Eigen::ArrayXpr>();
    ndarray::EigenView<double,2,-1,Eigen::ArrayXpr> out(output);
    if (!add) out.setZero();
    // We expect the Jacobian to be pretty sparse, so instead of just doing
    // standard multiplications here, we inspect each element and only do the
    // products we'll need.
    // This way if the user wants, for example, the derivatives with respect
    // to the ellipticity, we don't waste time computing elements that are
    // only useful when computing the derivatives wrt the centroid.
    double const eps = std::numeric_limits<double>::epsilon() * jacobian.lpNorm<Eigen::Infinity>();
    for (int n = 0; n < jacobian.cols(); ++n) {
        if (std::abs(jacobian(AT::XX, n)) > eps)
            out.col(n) += jacobian(AT::XX, n) * _xy.col(0).array() * dfdx;
        if (std::abs(jacobian(AT::XY, n)) > eps)
            out.col(n) += jacobian(AT::XY, n) * _xy.col(1).array() * dfdx;
        if (std::abs(jacobian(AT::X, n)) > eps)
            out.col(n) += jacobian(AT::X, n) * dfdx;
        if (std::abs(jacobian(AT::YX, n)) > eps)
            out.col(n) += jacobian(AT::YX, n) * _xy.col(0).array() * dfdy;
        if (std::abs(jacobian(AT::YY, n)) > eps)
            out.col(n) += jacobian(AT::YY, n) * _xy.col(1).array() * dfdy;
        if (std::abs(jacobian(AT::Y, n)) > eps)
            out.col(n) += jacobian(AT::Y, n) * dfdy;
    }
}

void GaussianModelBuilder::setOutput(ndarray::Array<double,1,1> const & array) {
    if (array.getSize<0>() != _xy.rows()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect size for array: got %d, expected %d")
             % array.getSize<0>() % _xy.rows()).str()
        );
    }
}

}}}} // namespace lsst::meas::extensions::multiShapelet
