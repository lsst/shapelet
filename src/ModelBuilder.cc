// -*- LSST-C++ -*-
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

#include "boost/format.hpp"
#include "boost/make_shared.hpp"

#include "lsst/utils/PowFast.h"
#include "lsst/pex/exceptions.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "ndarray/eigen.h"

namespace lsst { namespace shapelet {

namespace {

utils::PowFast const & powFast = utils::getPowFast<11>();

struct PowFastExpFunctor {
    inline float operator()(float x) const {
        return powFast.exp(x);
    }
};

template <typename T>
void fillExp(
    Eigen::Array<T,Eigen::Dynamic,1> const & x, Eigen::Array<T,Eigen::Dynamic,1> const & y,
    Eigen::Array<T,Eigen::Dynamic,1> & workspace,
    bool useApproximateExp
) {
    if (useApproximateExp) {
        workspace = (-0.5 * (x.square() + y.square())).unaryExpr(PowFastExpFunctor());
    } else {
        workspace = (-0.5 * (x.square() + y.square())).exp();
    }
}

template <typename T>
void fillHermite1d(
    int order,
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & workspace,
    Eigen::Array<T,Eigen::Dynamic,1> const & coord
) {
    if (order >= workspace.cols()) {
        workspace.resize(coord.size(), order + 1);
    }
    if (workspace.cols() > 0) {
        workspace.col(0).setConstant(BASIS_NORMALIZATION);
    }
    if (workspace.cols() > 1) {
        workspace.col(1) = intSqrt(2) * coord * workspace.col(0);
    }
    for (int j = 2; j <= order; ++j) {
        workspace.col(j) = rationalSqrt(2, j) * coord * workspace.col(j-1)
            - rationalSqrt(j - 1, j) * workspace.col(j-2);
    }
}

} // anonymous

template <typename T>
ModelBuilder<T>::ModelBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    bool useApproximateExp
) : _wsOrder(-1), _useApproximateExp(useApproximateExp),
    _ellipseFactor(1.0), _x(x), _y(y),
    _xt(_x.size()), _yt(_y.size()), _expWorkspace(_x.size())
{
    LSST_THROW_IF_NE(
        _x.size(), _y.size(),
        pex::exceptions::LengthErrorException,
        "x (%d) and y (%d) array sizes do not match"
    );
}

template <typename T>
void ModelBuilder<T>::update(afw::geom::ellipses::BaseCore const & ellipse) {
    afw::geom::ellipses::BaseCore::GridTransform gt(ellipse);
    typedef afw::geom::LinearTransform LT;
    LT transform = gt;
    _xt = _x * transform[LT::XX] + _y * transform[LT::XY];
    _yt = _x * transform[LT::YX] + _y * transform[LT::YY];
    _ellipseFactor = gt.getDeterminant();
    _wsOrder = -1;
    fillExp(_xt, _yt, _expWorkspace, _useApproximateExp);
}

template <typename T>
void ModelBuilder<T>::update(afw::geom::ellipses::Ellipse const & ellipse) {
    afw::geom::ellipses::Ellipse::GridTransform gt(ellipse);
    typedef afw::geom::AffineTransform AT;
    AT transform = gt;
    _xt = _x * transform[AT::XX] + _y * transform[AT::XY] + transform[AT::X];
    _yt = _x * transform[AT::YX] + _y * transform[AT::YY] + transform[AT::Y];
    _ellipseFactor = gt.getDeterminant();
    _wsOrder = -1;
    fillExp(_xt, _yt, _expWorkspace, _useApproximateExp);
}

template <typename T>
void ModelBuilder<T>::addModelMatrix(int order, ndarray::Array<T,2,-1> const & output) {
    if (output.template getSize<1>() != computeSize(order)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of columns of output matrix (%d) does not match shapelet order (%d->%d)")
             % output.template getSize<1>() % order % computeSize(order)).str()
        );
    }
    if (output.template getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of rows of output matrix (%d) does not match coordinate array (%d)")
             % output.template getSize<0>() % _x.size()).str()
        );
    }
    if (_wsOrder < order) {
        fillHermite1d(order, _xWorkspace, _xt);
        fillHermite1d(order, _yWorkspace, _yt);
        _wsOrder = order;
    }
    ndarray::EigenView<T,2,-1,Eigen::ArrayXpr> model(output);
    for (PackedIndex i; i.getOrder() <= order; ++i) {
        model.col(i.getIndex()) += _ellipseFactor * _expWorkspace
            * _xWorkspace.col(i.getX()) * _yWorkspace.col(i.getY());
    }
}

template <typename T>
void ModelBuilder<T>::addModelVector(
    int order,
    ndarray::Array<T const,1,1> const & coefficients,
    ndarray::Array<T,1,1> const & output
) {
    if (coefficients.template getSize<0>() != computeSize(order)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of coefficients (%d) does not match shapelet order (%d->%d)")
             % coefficients.template getSize<0>() % order % computeSize(order)).str()
        );
    }
    if (output.template getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of rows of output matrix (%d) does not match coordinate array (%d)")
             % output.template getSize<0>() % _x.size()).str()
        );
    }
    if (_wsOrder < order) {
        fillHermite1d(order, _xWorkspace, _xt);
        fillHermite1d(order, _yWorkspace, _yt);
        _wsOrder = order;
    }
    ndarray::EigenView<T,1,1,Eigen::ArrayXpr> model(output);
    for (PackedIndex i; i.getOrder() <= order; ++i) {
        model += _ellipseFactor * coefficients[i.getIndex()] * _expWorkspace
            * _xWorkspace.col(i.getX()) * _yWorkspace.col(i.getY());
    }
}

template class ModelBuilder<float>;
template class ModelBuilder<double>;

}} // namespace lsst::shapelet
