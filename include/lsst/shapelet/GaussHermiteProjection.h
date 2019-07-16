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

#ifndef LSST_AFW_MATH_SHAPELETS_GaussHermiteProjection
#define LSST_AFW_MATH_SHAPELETS_GaussHermiteProjection

#include "Eigen/Core"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/geom.h"
#include "lsst/shapelet/HermiteTransformMatrix.h"

namespace lsst { namespace shapelet {

class GaussHermiteProjection {
public:

    /// @brief Compute a matrix that projects from one shapelet basis ellipse to another.
    Eigen::MatrixXd compute(
        afw::geom::ellipses::Quadrupole const & inputEllipse, int inputOrder,
        afw::geom::ellipses::Quadrupole const & outputEllipse, int outputOrder
    ) const {
        return compute(
            inputEllipse.getGridTransform().getMatrix(), inputOrder,
            outputEllipse.getGridTransform().getMatrix(), outputOrder
        );
    }

    /// @brief Compute a matrix that projects from one shapelet basis "grid transform" to another.
    Eigen::MatrixXd compute(
        geom::LinearTransform const & inputTransform, int inputOrder,
        geom::LinearTransform const & outputTransform, int outputOrder
    ) const {
        return compute(
            inputTransform.getMatrix(), inputOrder,
            outputTransform.getMatrix(), outputOrder
        );
    }

    /// @brief Compute a matrix that projects from one shapelet basis "grid transform" to another.
    Eigen::MatrixXd compute(
        Eigen::Matrix2d const & inputTransform, int inputOrder,
        Eigen::Matrix2d const & outputTransform, int outputOrder
    ) const;

    int getMaxOrder() const { return _htm.getOrder(); }

    explicit GaussHermiteProjection(int maxOrder) : _htm(maxOrder) {}

private:
    HermiteTransformMatrix _htm;
};
}} // namespace lsst::shapelet

#endif // !LSST_AFW_MATH_SHAPELETS_GaussHermiteProjection
