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
#include "boost/make_shared.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "ndarray/eigen.h"

namespace lsst { namespace shapelet {

namespace {

void fillHermite1d(int order, Eigen::ArrayXXd & workspace, Eigen::ArrayXd const & coord) {
    if (order >= workspace.cols()) {
        workspace.resize(coord.size(), order + 1);
    }
    if (workspace.cols() > 0)
        workspace.col(0) = BASIS_NORMALIZATION * (-0.5 * coord.square()).exp();
    if (workspace.cols() > 1)
        workspace.col(1) = std::sqrt(2.0) * coord * workspace.col(0);
    for (int j = 2; j <= order; ++j) {
        workspace.col(j) = std::sqrt(2.0 / j) * coord * workspace.col(j-1)
            - std::sqrt((j - 1.0) / j) * workspace.col(j-2);
    }
}

} // anonymous

ModelBuilder::ModelBuilder(
    ndarray::Array<double const,1,1> const & x,
    ndarray::Array<double const,1,1> const & y
) : _wsOrder(-1),
    _x(x), _y(y),
    _xt(_x.size()), _yt(_y.size())
{}

void ModelBuilder::update(afw::geom::ellipses::BaseCore const & ellipse) {
    typedef afw::geom::LinearTransform LT;
    LT transform = ellipse.getGridTransform();
    _xt = _x * transform[LT::XX] + _y * transform[LT::XY];
    _yt = _x * transform[LT::YX] + _y * transform[LT::YY];
    _wsOrder = -1;
}

void ModelBuilder::addModelMatrix(int order, ndarray::Array<double,2,-1> const & output) {
    if (output.getSize<1>() != computeSize(order)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of columns of output matrix (%d) does not match shapelet order (%d->%d)")
             % output.getSize<1>() % order % computeSize(order)).str()
        );
    }
    if (output.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of rows of output matrix (%d) does not match coordinate array (%d)")
             % output.getSize<0>() % _x.size()).str()
        );
    }
    if (_wsOrder < order) {
        fillHermite1d(order, _xWorkspace, _xt);
        fillHermite1d(order, _yWorkspace, _yt);
        _wsOrder = order;
    }
    ndarray::EigenView<double,2,-1,Eigen::ArrayXpr> model(output);
    for (PackedIndex i; i.getOrder() <= order; ++i) {
        model.col(i.getIndex()) += _xWorkspace.col(i.getX()) * _yWorkspace.col(i.getY());
    }
}

void ModelBuilder::addModelVector(
    int order,
    ndarray::Array<double const,1,1> const & coefficients,
    ndarray::Array<double,1,1> const & output
) {
    if (coefficients.getSize<0>() != computeSize(order)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of coefficients (%d) does not match shapelet order (%d->%d)")
             % coefficients.getSize<0>() % order % computeSize(order)).str()
        );
    }
    if (output.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Number of rows of output matrix (%d) does not match coordinate array (%d)")
             % output.getSize<0>() % _x.size()).str()
        );
    }
    if (_wsOrder < order) {
        fillHermite1d(order, _xWorkspace, _xt);
        fillHermite1d(order, _yWorkspace, _yt);
        _wsOrder = order;
    }
    ndarray::EigenView<double,1,1,Eigen::ArrayXpr> model(output);
    for (PackedIndex i; i.getOrder() <= order; ++i) {
        model += coefficients[i.getIndex()] * _xWorkspace.col(i.getX()) * _yWorkspace.col(i.getY());
    }
}

}} // namespace lsst::shapelet
