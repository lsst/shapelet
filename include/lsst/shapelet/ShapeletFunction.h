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
 
#ifndef LSST_AFW_MATH_SHAPELETS_SHAPELETFUNCTION_H
#define LSST_AFW_MATH_SHAPELETS_SHAPELETFUNCTION_H

#include "ndarray.h"
#include "lsst/shapelet/constants.h"
#include "lsst/shapelet/GaussHermiteEvaluator.h"
#include "lsst/shapelet/ConversionMatrix.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/image/Image.h"

#include <list>

namespace lsst { namespace shapelet {

class ShapeletFunctionEvaluator;

/**
 *  @brief A 2-d function defined by an expansion onto a Gauss-Laguerre or Gauss-Hermite basis.
 *
 *  The coefficients are in units of flux, not surface brightness; increasing
 *  the area of the basis ellipse while leaving the coefficients unchanged will
 *  decrease the surface brightness by the ratio of the areas of the new and old
 *  ellipses.  This convention is necessary to ensure that convolution with a
 *  delta function is sane.  Because the BasisEvaluator class does not deal with
 *  ellipses, it is necessary to divide its output by R^2 to get the same result
 *  as a ShapeletFunctionEvaluator whose ellipse has determinant radius of R.
 *
 *  The coefficient normalization is not identical to that of a Gaussian, however;
 *  a zeroth-order ShapeletFunction with its only coefficient value set to 1 has a
 *  flux of 2.0 * pi^(1/2).  This value is defined as ShapeletFunction::FLUX_FACTOR.
 *  Of course, to get the flux of a more complex shapelet expansion you have to use
 *  ShapeletFunctionEvaluator::integrate().
 *
 *  Note that the units of the coefficients would have to be radius for basis
 *  functions with the same ellipse to be orthonormal, but this orthonormality
 *  isn't very useful, because the basis functions aren't even orthogonal in the
 *  more common case that the ellipses differ.  However, basis functions defined
 *  on the unit circle are still orthonormal.
 */
class ShapeletFunction {
public:

    typedef boost::shared_ptr<ShapeletFunction> Ptr;
    typedef boost::shared_ptr<ShapeletFunction const> ConstPtr;

    typedef ShapeletFunctionEvaluator Evaluator;

    static double const FLUX_FACTOR;

    /// @brief Return the maximum order (inclusive), either @f$n_x + n_y@f$ or @f$p + q@f$.
    int getOrder() const { return _order; }

    /// @brief Get the ellipse (const).
    afw::geom::ellipses::Ellipse const & getEllipse() const { return _ellipse; }

    /// @brief Get the ellipse (non-const).
    afw::geom::ellipses::Ellipse & getEllipse() { return _ellipse; }

    /// @brief Set the ellipse.
    void setEllipse(afw::geom::ellipses::Ellipse const & ellipse) { _ellipse = ellipse; }
    
    /// @brief Return the basis type (HERMITE or LAGUERRE).
    BasisTypeEnum getBasisType() const { return _basisType; }

    /// @brief Change the basis type and convert coefficients in-place correspondingly.
    void changeBasisType(BasisTypeEnum basisType) {
        ConversionMatrix::convertCoefficientVector(_coefficients, _basisType, basisType, _order);
        _basisType = basisType;
    }

    /// @brief Normalize the integral of the shapelet function to 1.
    void normalize();

    /// @brief Return the coefficient vector.
    ndarray::Array<double,1,1> const getCoefficients() { return _coefficients; }

    /// @brief Return the coefficient vector (const).
    ndarray::Array<double const,1,1> const getCoefficients() const { return _coefficients; }

    /// @brief Convolve the shapelet function.
    ShapeletFunction convolve(ShapeletFunction const & other) const;

    /// @brief Construct a helper object that can efficiently evaluate the function.
    Evaluator evaluate() const;

    /// @brief Shift the shapelet function by shifting the basis ellipse.
    void shiftInPlace(afw::geom::Extent2D const & offset) {
        _ellipse.getCenter() += offset;
    }

    /// @brief Transform the shapelet function by transforming the basis ellipse.
    void transformInPlace(afw::geom::AffineTransform const & transform) {
        _ellipse.transform(transform).inPlace();
    }

    /// @brief Construct a function with a unit-circle ellipse and set all coefficients to zero.
    ShapeletFunction(int order, BasisTypeEnum basisType);

    /// @brief Construct a function with a unit-circle ellipse and a deep-copied coefficient vector.
    ShapeletFunction(
        int order, BasisTypeEnum basisType,
        ndarray::Array<double,1,1> const & coefficients
    );

    /// @brief Construct a function with a circular ellipse and set all coefficients to zero.
    ShapeletFunction(int order, BasisTypeEnum basisType, double radius,
                     afw::geom::Point2D const & center=afw::geom::Point2D());

    /// @brief Construct a function with a circular ellipse and a deep-copied coefficient vector.
    ShapeletFunction(
        int order, BasisTypeEnum basisType, double radius, afw::geom::Point2D const & center,
        ndarray::Array<double,1,1> const & coefficients
    );

    /// @brief Construct a function and set all coefficients to zero.
    ShapeletFunction(int order, BasisTypeEnum basisType,
        afw::geom::ellipses::Ellipse const & ellipse);

    /// @brief Construct a function with a deep-copied coefficient vector.
    ShapeletFunction(
        int order, BasisTypeEnum basisType,
        afw::geom::ellipses::Ellipse const & ellipse,
        ndarray::Array<double const,1,1> const & coefficients
    );

    /// @brief Copy constructor (deep).
    ShapeletFunction(ShapeletFunction const & other);

    /// @brief Assignment (deep).
    ShapeletFunction & operator=(ShapeletFunction const & other);

    /// @brief Default constructor to appease SWIG (used by std::list).  Not for use by users.
    ShapeletFunction();

private:

    int _order;
    BasisTypeEnum _basisType;
    lsst::afw::geom::ellipses::Ellipse _ellipse;
    ndarray::Array<double,1,1> _coefficients;
};

/**
 *  @brief Evaluates a ShapeletFunction.
 *
 *  This is distinct from ShapeletFunction to allow the evaluator to construct temporaries
 *  and allocate workspace that will be reused when evaluating at multiple points.
 *
 *  A ShapeletFunctionEvaluator is invalidated whenever the ShapeletFunction it
 *  was constructed from is modified.
 */
class ShapeletFunctionEvaluator {
public:

    typedef boost::shared_ptr<ShapeletFunctionEvaluator> Ptr;
    typedef boost::shared_ptr<ShapeletFunctionEvaluator const> ConstPtr;

    /// @brief Evaluate at the given point.
    double operator()(double x, double y) const {
        return this->operator()(afw::geom::Point2D(x, y));
    }

    /// @brief Evaluate at the given point.
    double operator()(afw::geom::Point2D const & point) const {
        return _normalization * _h.sumEvaluation(_coefficients, _transform(point));
    }

    /// @brief Evaluate at the given point.
    double operator()(afw::geom::Extent2D const & point) const {
        return _normalization * _h.sumEvaluation(_coefficients, _transform(point));
    }

    /// @brief Evaluate at the given points, returning a newly-allocated array.
    ndarray::Array<double,1,1> operator()(
        ndarray::Array<double const,1> const & x,
        ndarray::Array<double const,1> const & y
    ) const;

    /// @brief Add the function to the given image-like array.
    void addToImage(
        ndarray::Array<double,2,1> const & array,
        afw::geom::Point2I const & xy0 = afw::geom::Point2I()
    ) const;

    /// @brief Evaluate the function on the given image.
    void addToImage(afw::image::Image<double> & image) const {
        addToImage(image.getArray(), image.getXY0());
    }

    /// @brief Compute the definite integral or integral moments.
    double integrate() const {
        return _h.sumIntegration(_coefficients);
    }

    /// @brief Return the unweighted dipole and quadrupole moments of the function as an ellipse.
    afw::geom::ellipses::Ellipse computeMoments() const;

    /// @brief Update the evaluator from the given function.
    void update(ShapeletFunction const & function);

    /// @brief Construct an evaluator for the given function.
    explicit ShapeletFunctionEvaluator(ShapeletFunction const & function);

private:
    
    friend class MultiShapeletFunctionEvaluator;

    void _computeRawMoments(double & q0, Eigen::Vector2d & q1, Eigen::Matrix2d & q2) const;

    double _normalization;
    ndarray::Array<double const,1,1> _coefficients;
    afw::geom::AffineTransform _transform;
    GaussHermiteEvaluator _h;
};

inline ShapeletFunctionEvaluator ShapeletFunction::evaluate() const {
    return ShapeletFunctionEvaluator(*this);
}

}} // namespace lsst::shapelet

#endif // !defined(LSST_AFW_MATH_SHAPELETS_SHAPELETFUNCTION_H)
