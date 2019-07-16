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

#ifndef LSST_AFW_MATH_SHAPELETS_HERMITEEVALUATOR_H
#define LSST_AFW_MATH_SHAPELETS_HERMITEEVALUATOR_H

#include "ndarray.h"
#include "lsst/afw/geom.h"
#include "lsst/shapelet/constants.h"
#include "Eigen/Core"

namespace lsst { namespace shapelet {

/**
 *  @brief An iterator-like object to help in traversing "packed" shapelet or Hermite polynomial
 *         matrix or vector dimensions.
 *
 *  A pair of indices (x,y) is mapped to the packed position i = (x+y)(x+y+1)/2 + x.
 *
 *  Typical usage is in a nested loop of the form:
 *  @code
 *      for (PackedIndex i; i.getOrder() <= order; ++i) {
 *          // utilize i
 *      }
 *  @endcode
 */
class PackedIndex {
public:

    static int const computeOffset(int order) { return order*(order+1)/2; }
    static int const computeIndex(int x, int y) { return computeOffset(x+y) + x; }

    PackedIndex & operator++() {
        ++_i;
        if (--_y < 0) {
            _x = 0;
            _y = ++_n;
        } else {
            ++_x;
        }
        return *this;
    }

    int const getOrder() const { return _n; }
    int const getX() const { return _x; }
    int const getY() const { return _y; }

    int const getIndex() const { return _i; }

    PackedIndex() : _n(0), _i(0), _x(0), _y(0) {}
    PackedIndex(int const x, int const y) : _n(x+y), _i(computeOffset(_n) + x), _x(x), _y(y) {}

private:
    int _n;
    int _i;
    int _x;
    int _y;
};

/**
 *  @brief A class to evaluate HERMITE shapelet-related quantities.
 */
class GaussHermiteEvaluator {
public:

    /**
     *  @brief Fill a matrix with the function inner products of two
     *         HERMITE shapelet basis functions with different scales.
     *  @f$
     *      M_{\mathbf{i},\mathbf{j}} =
     *      \int d^2 \mathbf{x} \psi_\mathbf{i}(a\mathbf{x})\phi_\mathbf{j}(b\mathbf{x})
     *  @f$
     */
    static Eigen::MatrixXd computeInnerProductMatrix(
        int rowOrder, int colOrder, double a, double b
    );

    int getOrder() const { return _xWorkspace.getSize<0>() - 1; }

    /**
     *  @brief Fill a vector whose dot product with a HERMITE coefficient vector evaluates a
     *         simple unscaled shapelet expansion at the given point.
     */
    void fillEvaluation(
        Array1d const & target, double x, double y,
        Array1d const & dx = Array1d(),
        Array1d const & dy = Array1d()
    ) const;

    /**
     *  @brief Fill a vector whose dot product with a HERMITE coefficient vector evaluates a
     *         simple unscaled shapelet expansion at the given point.
     */
    void fillEvaluation(
        Array1d const & target, afw::geom::Point2D const & point,
        Array1d const & dx = Array1d(),
        Array1d const & dy = Array1d()
    ) const {
        fillEvaluation(target, point.getX(), point.getY(), dx, dy);
    }

    /**
     *  @brief Fill a vector whose dot product with a HERMITE coefficient vector evaluates a
     *         simple unscaled shapelet expansion at the given point.
     */
    void fillEvaluation(
        Array1d const & target, afw::geom::Extent2D const & point,
        Array1d const & dx = Array1d(),
        Array1d const & dy = Array1d()
    ) const {
        fillEvaluation(target, point.getX(), point.getY(), dx, dy);
    }

    /**
     *  @brief Fill a vector whose dot product with a HERMITE coefficient vector integrates
     *         a simple unscaled shapelet expansion.
     */
    void fillIntegration(Array1d const & target, int xMoment=0, int yMoment=0) const;

    /**
     *  @brief Evaluate a simple unscaled shapelet expansion at the given point.
     */
    double sumEvaluation(
        ndarray::Array<double const,1> const & coeff, double x, double y,
        double * dx = 0, double * dy = 0
    ) const;

    /**
     *  @brief Evaluate a simple unscaled shapelet expansion at the given point.
     */
    double sumEvaluation(
        ndarray::Array<double const,1> const & coeff, afw::geom::Point2D const & point,
        double * dx = 0, double * dy = 0
    ) const {
        return sumEvaluation(coeff, point.getX(), point.getY(), dx, dy);
    }

    /**
     *  @brief Evaluate a simple unscaled shapelet expansion at the given point.
     */
    double sumEvaluation(
        ndarray::Array<double const,1> const & coeff, afw::geom::Extent2D const & point,
        double * dx = 0, double * dy = 0
    ) const {
        return sumEvaluation(coeff, point.getX(), point.getY(), dx, dy);
    }

    /**
     *  @brief Integrate a simple unscaled shapelet expansion.
     */
    double sumIntegration(ndarray::Array<double const,1> const & coeff, int xMoment=0, int yMoment=0) const;

    explicit GaussHermiteEvaluator(int order);

private:

    ndarray::Array<double,1,1> _xWorkspace;
    ndarray::Array<double,1,1> _yWorkspace;
    ndarray::Array<double,1,1> _dxWorkspace;
    ndarray::Array<double,1,1> _dyWorkspace;
};

}} // namespace lsst::shapelet

#endif // !defined(LSST_AFW_MATH_SHAPELETS_HERMITEEVALUATOR_H)
