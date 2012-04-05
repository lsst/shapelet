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
#ifndef MULTISHAPELET_GaussianModelBuilder_h_INCLUDED
#define MULTISHAPELET_GaussianModelBuilder_h_INCLUDED

#include "ndarray/eigen.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class GaussianModelBuilder {
public:

    explicit GaussianModelBuilder(afw::detection::Footprint const & region);

    explicit GaussianModelBuilder(afw::geom::Box2I const & region);
    
    void update(afw::geom::ellipses::Ellipse const & ellipse);
    
    ndarray::Array<double const,1,1> getModel() const { return _model; }

    int getSize() const { return _xy.rows(); }

    void computeDerivative(
        ndarray::Array<double,2,-1> const & output
    ) const;

    void computeDerivative(
        ndarray::Array<double,2,-1> const & output,
        Eigen::Matrix<double,5,Eigen::Dynamic> const & jacobian,
        bool add = false 
    ) const;

    void setOutput(ndarray::Array<double,1,1> const & array);

private:

    void _computeDerivative(
        ndarray::Array<double,2,-1> const & output,
        Eigen::Matrix<double,6,Eigen::Dynamic> const & jacobian,
        bool add
    ) const;

    PTR(afw::geom::ellipses::Ellipse) _ellipse;
    ndarray::Array<double,1,1> _model;
    Eigen::Matrix<double,Eigen::Dynamic,2> _xy;
    Eigen::Matrix<double,Eigen::Dynamic,2> _xyt;
};


}}}} // namespace lsst::meas::extensions::multisShapelet

#endif // !MULTISHAPELET_GaussianModelBuilder_h_INCLUDED
