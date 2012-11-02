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

#include "lsst/shapelet/GaussHermiteProjection.h"

namespace lsst { namespace shapelet {

Eigen::MatrixXd GaussHermiteProjection::compute(
    Eigen::Matrix2d const & inputTransform, int inputOrder,
    Eigen::Matrix2d const & outputTransform, int outputOrder
) const {
    assert(inputOrder <= getMaxOrder());
    assert(outputOrder <= getMaxOrder());
    int fullOrder = std::max(inputOrder, outputOrder);
    int inputSize = computeSize(inputOrder);
    int outputSize = computeSize(outputOrder);
    Eigen::Matrix2d wwtInv = 0.5 * (
        inputTransform.adjoint() * inputTransform
        + outputTransform.adjoint() * outputTransform
    );
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(wwtInv);
    Eigen::Matrix2d w = eig.operatorInverseSqrt();
    Eigen::MatrixXd q1 = _htm.compute(outputTransform * w, fullOrder);
    Eigen::MatrixXd q2 = _htm.compute(inputTransform * w, fullOrder);
    Eigen::MatrixXd result = q1.adjoint().topRows(outputSize) * q2.leftCols(inputSize);
    result *= outputTransform.determinant() / std::sqrt(eig.eigenvalues().prod());
    return result;
}

}} // namespace lsst::shapelet
