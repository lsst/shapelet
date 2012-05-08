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

void fillHermite1d(Eigen::ArrayXXd & workspace, Eigen::ArrayXd const & coord) {
    if (workspace.cols() > 0)
        workspace.col(0) = NORMALIZATION * (-0.5 * coord.square()).exp();
    if (workspace.cols() > 1)
        workspace.col(1) = std::sqrt(2.0) * coord * workspace.col(0);
    for (int j = 2; j < workspace.cols(); ++j) {
        workspace.col(j) = std::sqrt(2.0 / j) * coord * workspace.col(j-1)
            - std::sqrt((j - 1.0) / j) * workspace.col(j-2);
    }
}

} // anonymous

ModelBuilder::ModelBuilder(
    int order,
    ndarray::Array<double const,1,1> const & x,
    ndarray::Array<double const,1,1> const & y
) : _order(order),
    _model(ndarray::allocate(x.getSize<0>(), computeSize(order))),
    _x(x), _y(y),
    _xt(_x.size()), _yt(_y.size()),
    _xWorkspace(_x.size(), order + 1), _yWorkspace(_x.size(), order + 1)
{}

void ModelBuilder::update(afw::geom::ellipses::BaseCore const & ellipse) {
    typedef afw::geom::LinearTransform LT;
    LT transform = ellipse.getGridTransform();
    _xt = _x * transform[LT::XX] + _y * transform[LT::XY];
    _yt = _x * transform[LT::YX] + _y * transform[LT::YY];
    fillHermite1d(_xWorkspace, _xt);
    fillHermite1d(_yWorkspace, _yt);
    ndarray::EigenView<Pixel,2,-2,Eigen::ArrayXpr> model(_model);
    for (PackedIndex i; i.getOrder() <= _order; ++i) {
        model.col(i.getIndex()) = _xWorkspace.col(i.getX()) * _yWorkspace.col(i.getY());
    }
}

}} // namespace lsst::shapelet
