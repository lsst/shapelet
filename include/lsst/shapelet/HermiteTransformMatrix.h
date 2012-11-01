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

#ifndef LSST_AFW_MATH_SHAPELETS_HermiteTransformMatrix
#define LSST_AFW_MATH_SHAPELETS_HermiteTransformMatrix

#include "Eigen/Core"

#include "lsst/afw/geom.h"
#include "lsst/shapelet/constants.h"

namespace lsst { namespace shapelet {

/**
 *  @brief A class that computes a matrix that applies a linear transform to a 2-d Hermite
 *         polynomial expansion.
 *
 *  Let
 *  @f[
 *      Z_{\bm{n}}\!(\bm{x}) \equiv \mathcal{H}_{n_0}\!(x_0)\;\mathcal{H}_{n_1}\!(x_1)
 *  @f]
 *  where
 *  @f[
 *      \mathcal{H}_n(x)=(2^n \pi^{1/2} n!)^{-1/2}H_n(x)
 *  @f]
 *  is the @f$i@f$th "alternate" Hermite
 *  polynomial.  This function computes the matrix @f$\bm{Q}(\bm{U})@f$ given a linear
 *  transform @f$\bm{U}@f$ such that
 *  @f[
 *      Z_{\bm{m}}\!(\bm{U}\bm{x}) = \sum_{\bm{n}}Q_{\bm{m},\bm{n}}\!(\bm{U})\,Z_{\bm{n}}\!(\bm{x})
 *  @f]
 */
class HermiteTransformMatrix {
public:

    /// @brief Compute the matrix for a new linear transform.
    Eigen::MatrixXd compute(Eigen::Matrix2d const & transform) const {
        return compute(transform, _order);
    }

    /// @brief Compute the matrix for a new linear transform.
    Eigen::MatrixXd compute(afw::geom::LinearTransform const & transform) const {
        return compute(transform.getMatrix(), _order);
    }

    /// @brief Compute the matrix for a new linear transform at the given order (must be <= getOrder()).
    Eigen::MatrixXd compute(Eigen::Matrix2d const & transform, int order) const;

    /// @brief Compute the matrix for a new linear transform at the given order (must be <= getOrder()).
    Eigen::MatrixXd compute(afw::geom::LinearTransform const & transform, int order) const {
        return compute(transform.getMatrix(), order);
    }

    /**
     *  @brief Return the matrix that maps (1-d) regular polynomials to the alternate Hermite polynomials.
     *
     *  The matrix is always lower triangular, and has size equal to getOrder()+1.
     */
    Eigen::MatrixXd getCoefficientMatrix() const { return _coeffFwd; }

    /**
     *  @brief Return the matrix that maps (1-d) alternate Hermite polynomials to regular polynomials.
     *
     *  The matrix is always lower triangular, and has size equal to getOrder()+1.
     */
    Eigen::MatrixXd getInverseCoefficientMatrix() const { return _coeffInv; }

    /// @brief Return the maximum order at which the matrix can be computed.
    int getOrder() const { return _order; }

    /// @brief Construct an instance able to compute the transform matrix at up to the given order.
    explicit HermiteTransformMatrix(int order);

private:
    int _order;
    Eigen::MatrixXd _coeffFwd;
    Eigen::MatrixXd _coeffInv;
};

}} // namespace lsst::shapelet

#endif // !LSST_AFW_MATH_SHAPELETS_HermiteTransformMatrix
