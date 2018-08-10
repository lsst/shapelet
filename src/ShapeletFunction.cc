// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#include "boost/format.hpp"

#include "ndarray/eigen.h"
#include "lsst/shapelet/ShapeletFunction.h"
#include "lsst/shapelet/ConversionMatrix.h"
#include "lsst/shapelet/GaussHermiteConvolution.h"
#include "lsst/pex/exceptions.h"

namespace lsst { namespace shapelet {

namespace {

static inline void validateSize(int expected, int actual) {
    if (expected != actual) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format(
                "Coefficient vector for ShapeletFunction has incorrect size (%n, should be %n)."
            ) % actual % expected).str()
        );
    }
}

} // anonymous

double const ShapeletFunction::FLUX_FACTOR = 2.0 * std::sqrt(afw::geom::PI);

ShapeletFunction::ShapeletFunction() : 
    _order(0), _basisType(HERMITE),
    _ellipse(EllipseCore(1.0, 1.0, 0.0), afw::geom::Point2D()), 
    _coefficients(ndarray::allocate(1))
{
    _coefficients[0] = 0.0;
}

ShapeletFunction::ShapeletFunction(int order, BasisTypeEnum basisType) :
    _order(order), _basisType(basisType), _ellipse(EllipseCore()),
    _coefficients(ndarray::allocate(computeSize(_order)))
{
    _coefficients.deep() = 0.0;
}

ShapeletFunction::ShapeletFunction(
    int order, BasisTypeEnum basisType,
    ndarray::Array<double,1,1> const & coefficients
) :
    _order(order), _basisType(basisType), _ellipse(EllipseCore()),
    _coefficients(ndarray::copy(coefficients))
{
    validateSize(computeSize(order), _coefficients.getSize<0>());
}
 
ShapeletFunction::ShapeletFunction(
    int order, BasisTypeEnum basisType, double radius,
    afw::geom::Point2D const & center
) :
    _order(order), _basisType(basisType), _ellipse(EllipseCore(), center),
    _coefficients(ndarray::allocate(computeSize(_order)))
{
    _ellipse.getCore().scale(radius);
    _coefficients.deep() = 0.0;
}

ShapeletFunction::ShapeletFunction(
    int order, BasisTypeEnum basisType, double radius, afw::geom::Point2D const & center,
    ndarray::Array<double,1,1> const & coefficients
) :
    _order(order), _basisType(basisType), _ellipse(EllipseCore(), center),
    _coefficients(ndarray::copy(coefficients))
{
    _ellipse.getCore().scale(radius);
    validateSize(computeSize(order), _coefficients.getSize<0>());
}
 
ShapeletFunction::ShapeletFunction(
    int order, BasisTypeEnum basisType, afw::geom::ellipses::Ellipse const & ellipse
) :
    _order(order), _basisType(basisType), _ellipse(EllipseCore(ellipse.getCore()), ellipse.getCenter()),
    _coefficients(ndarray::allocate(computeSize(_order)))
{
    _coefficients.deep() = 0.0;
}

ShapeletFunction::ShapeletFunction(
    int order, BasisTypeEnum basisType, afw::geom::ellipses::Ellipse const & ellipse,
    ndarray::Array<double const,1,1> const & coefficients
) :
    _order(order), _basisType(basisType), _ellipse(EllipseCore(ellipse.getCore()), ellipse.getCenter()),
    _coefficients(ndarray::copy(coefficients))
{
    validateSize(computeSize(order), _coefficients.getSize<0>());
}

ShapeletFunction::ShapeletFunction(ShapeletFunction const & other) :
    _order(other._order), _basisType(other._basisType), _ellipse(other._ellipse),
    _coefficients(ndarray::copy(other._coefficients))
{}

ShapeletFunction & ShapeletFunction::operator=(ShapeletFunction const & other) {
    if (&other != this) {
        if (other.getOrder() != this->getOrder()) {
            _order = other.getOrder();
            _coefficients = ndarray::copy(other.getCoefficients());
        } else {
            _coefficients.deep() = other.getCoefficients();
        }
        _basisType = other.getBasisType();
        _ellipse = other.getEllipse();
    }
    return *this;
}

void ShapeletFunction::normalize(double value) {
    _coefficients.deep() *= value / evaluate().integrate();
}

void ShapeletFunctionEvaluator::update(ShapeletFunction const & function) {
    validateSize(_h.getOrder(), function.getOrder());
    afw::geom::ellipses::Ellipse::GridTransform gt(function.getEllipse());
    _transform = gt;
    _normalization = gt.getDeterminant();
    switch (function.getBasisType()) {
    case HERMITE:
        _coefficients = function.getCoefficients();
        break;
    case LAGUERRE:
        ndarray::Array<double,1,1> tmp(ndarray::copy(function.getCoefficients()));
        ConversionMatrix::convertCoefficientVector(
            tmp, LAGUERRE, HERMITE, function.getOrder()
        );
        _coefficients = tmp;
        break;
    }
}

ShapeletFunctionEvaluator::ShapeletFunctionEvaluator(
    ShapeletFunction const & function
) : _normalization(1.0),
    _transform(),
    _h(function.getOrder())
{
    update(function);
}

ShapeletFunction ShapeletFunction::convolve(ShapeletFunction const & other) const {
    GaussHermiteConvolution convolution(getOrder(), other);
    afw::geom::ellipses::Ellipse newEllipse(_ellipse);
    auto matrixNdArray = convolution.evaluate(newEllipse);
    if (_basisType == LAGUERRE) {
        ConversionMatrix::convertCoefficientVector(_coefficients, LAGUERRE, HERMITE, getOrder());
    }
    ShapeletFunction result(convolution.getRowOrder(), HERMITE);
    result.setEllipse(newEllipse);
    ndarray::asEigenMatrix(result.getCoefficients()) =
            ndarray::asEigenMatrix(matrixNdArray) * ndarray::asEigenMatrix(_coefficients);
    if (_basisType == LAGUERRE) {
        ConversionMatrix::convertCoefficientVector(_coefficients, HERMITE, LAGUERRE, getOrder());
        result.changeBasisType(LAGUERRE);
    }
    return result;
}

ndarray::Array<double,1,1> ShapeletFunctionEvaluator::operator()(
    ndarray::Array<double const,1> const & x,
    ndarray::Array<double const,1> const & y
) const {
    ndarray::Array<double,1,1> output = ndarray::allocate(x.getSize<0>());
    for (int i = 0, n = x.getSize<0>(); i < n; ++i) {
        output[i] = (*this)(x[i], y[i]);
    }
    return output;
}

void ShapeletFunctionEvaluator::addToImage(
    ndarray::Array<double,2,1> const & array,
    afw::geom::Point2I const & xy0
) const {
    ndarray::Array<double,2,1>::Iterator yIter = array.begin();
    for (int y = xy0.getY(); yIter != array.end(); ++y, ++yIter) {
        ndarray::Array<double,2,1>::Reference::Iterator xIter = yIter->begin();
        for (int x = xy0.getX(); xIter != yIter->end(); ++x, ++xIter) {
            *xIter += (*this)(x, y);
        }
    }
}

void ShapeletFunctionEvaluator::_computeRawMoments(
    double & q0, Eigen::Vector2d & q1, Eigen::Matrix2d & q2
) const {
    Eigen::Matrix2d a = _transform.getLinear().inverted().getMatrix();
    Eigen::Vector2d b = _transform.getTranslation().asEigen();

    double m0 = _h.sumIntegration(_coefficients, 0, 0);
    q0 += m0;

    Eigen::Vector2d m1(
        _h.sumIntegration(_coefficients, 1, 0),
        _h.sumIntegration(_coefficients, 0, 1)
    );
    q1 += a * (m1 - b * m0);

    Eigen::Matrix2d m2;
    m2(0, 0) = _h.sumIntegration(_coefficients, 2, 0);
    m2(1, 1) = _h.sumIntegration(_coefficients, 0, 2);
    m2(0, 1) = m2(1, 0) = _h.sumIntegration(_coefficients, 1, 1);
    q2 += a * (m2 + b * b.transpose() * m0 - m1 * b.transpose() - b * m1.transpose()) * a.transpose();
}

afw::geom::ellipses::Ellipse ShapeletFunctionEvaluator::computeMoments() const {
    double q0 = 0.0;
    Eigen::Vector2d q1 = Eigen::Vector2d::Zero();
    Eigen::Matrix2d q2 = Eigen::Matrix2d::Zero();
    _computeRawMoments(q0, q1, q2);
    q1 /= q0;
    q2 /= q0;
    q2 -= q1 * q1.transpose();
    return afw::geom::ellipses::Ellipse(
        afw::geom::ellipses::Quadrupole(afw::geom::ellipses::Quadrupole::Matrix(q2), false),
        afw::geom::Point2D(q1)
    );
}

}} // namespace lsst::shapelet
