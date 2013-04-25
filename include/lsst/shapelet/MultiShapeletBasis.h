// -*- lsst-c++ -*-
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
#ifndef LSST_SHAPELET_MultiShapeletBasis_h_INCLUDED
#define LSST_SHAPELET_MultiShapeletBasis_h_INCLUDED

#include "boost/scoped_ptr.hpp"
#include "boost/noncopyable.hpp"

#include "lsst/shapelet/MultiShapeletFunction.h"

namespace lsst { namespace shapelet {

/**
 *  @brief Simple struct that represents one shapelet expansion in a MultiShapeletBasis
 *
 *  A MultiShapeletBasis is formed by the linear combination of several shapelet bases
 *  with different radii and common ellipticity; this represents a single shapelet basis
 *  within the MultiShapeletBasis.
 *
 *  @note This really ought to be an inner class, and should generally be referred to via
 *        the MultiShapeletBasis::Component typedef, but Swig doesn't handle inner classes.
 */
class MultiShapeletBasisComponent {
public:

    /**
     *  @brief Main constructor for MultiShapeletBasisComponent
     *
     *  Should usually only be called by MultiShapeletBasis::addComponent.
     *
     *  @param[in] radius   Radius of the shapelet expansion defined by this component.
     *  @param[in] order    Order of the shapelet expansion.
     *  @param[in] matrix   Matrix whose elements [i,j] map MultiShapeletBasis elements j to shapelet
     *                      terms i; must have dimensions [computeSize(order), basis.getSize()], where
     *                      "basis" is the MultiShapeletBasis this component is attached to.  Will be
     *                      deep-copied by the constructor.
     *
     *  Note that matrix elements follow the amplitude convention defined by ShapeletFunction; values
     *  are proportional to flux, not surface brightness.
     */
    MultiShapeletBasisComponent(double radius, int order, ndarray::Array<double const,2,2> const & matrix);

    /// Return the radius of this shapelet expansion
    double getRadius() const { return _radius; }

    /// Order of this shapelet expansion
    int getOrder() const { return _order; }

    /// Matrix whose elements [i,j] map MultiShapeletBasis elements j to shapelet terms i
    ndarray::Array<double const,2,2> getMatrix() const { return _matrix; }

private:
    friend class MultiShapeletBasis;

    double _radius;
    int _order;
    ndarray::Array<double const,2,2> _matrix;
};

/**
 *  @brief A basis formed from a linear combination of shapelet bases that differ only in radius.
 *
 *  A MultiShapeletBasis can have many "components" (shapelet bases with different orders and radii),
 *  which are mapped via matrices into one or more "elements".  It's common for a basis to have only
 *  one or two elements, representing a galaxy model that is a linear combination of one or two
 *  radial profiles.  It's also common for most components to be zeroth order (Gaussians), as higher-
 *  order shapelet terms don't provide much of an advantage when approximating axisymmetric functions.
 *
 *  MultiShapeletBasis itself provides the interface to define these multi-Gaussian approximations
 *  and combine and refine them, and delegates the work of defining them to MultiShapeletFunction
 *  (via the makeFunction() method) and the MultiShapeletMatrixBuilder class.  MultiShapeletFunction
 *  is a more user-friendly and versatile approach, intended for debugging and testing, while the
 *  MultiShapletMatrixBuilder approach is intended for performance-critical evaluation of PSF-convolved
 *  MultiShapeletBasis objects.
 */
class MultiShapeletBasis {
public:
    typedef MultiShapeletBasisComponent Component;
    typedef std::vector<Component> ComponentVector;
    typedef ComponentVector::const_iterator Iterator;

    /// Construct a MultiShapeletBasis with the given number of elements (i.e. free amplitudes).
    explicit MultiShapeletBasis(int size);

    /// Return the number of elements in the MultiShapeletBasis
    int getSize() const { return _size; }

    //@{
    /// Iterator over the components (distinct shapelet bases) of the MultiShapeletBasis
    Iterator begin() const { return _components.begin(); }
    Iterator end() const { return _components.end(); }
    //@}

    /**
     *  @brief Add a new component (shapelet basis) to the MultiShapeletBasis
     *
     *  @copydetails MultiShapeletBasisComponent::MultiShapeletBasisComponent
     */
    void addComponent(double radius, int order, ndarray::Array<double const,2,2> const & matrix);

    /// Multiply the radius of all basis elements by the given factor.
    void scale(double factor);

    /// Rescale all matrices so each element has unit flux.
    void normalize();

    /// Combine the given basis with this (in place), by appending its elements.
    void merge(MultiShapeletBasis const & other);

    /**
     *  @brief Create a MultiShapeletFunction from the basis.
     *
     *  @param[in]  ellipse       Shapelet basis ellipse that will define the MultiShapeletFunction
     *                            (will be scaled by the radius of each component).
     *  @param[in]  coefficients  Coefficients of the basis elements.
     */
    MultiShapeletFunction makeFunction(
        afw::geom::ellipses::Ellipse const & ellipse,
        ndarray::Array<double const,1,1> const & coefficients
    ) const;

private:
    int _size;
    ComponentVector _components;
};

/**
 *  @brief Class that efficiently evaluates a PSF-convolved elliptical MultiShapeletBasis matrix
 *         at predefined points.
 *
 *  The output "matrix" has pixels along the rows, and MultiShapeletBasis elements along columns;
 *  this is the design matrix involved in a linear least squares fit for the basis coefficients.
 */
template <typename T>
class MultiShapeletMatrixBuilder : private boost::noncopyable {
public:

    /**
     *  @brief Main constructor.
     *
     *  @param[in] basis     MultiShapeletBasis that defines the matrix to be generated
     *  @param[in] psf       Point-spread function with which to convolve the MultiShapeletBasis model
     *  @param[in] x         Array of x coordinates at which the matrix will be evaluated
     *  @param[in] y         Array of y coordinates at which the matrix will be evaluated (must have the
     *                       same size as x).
     *  @param[in] useApproximateExp   Use utils::PowFast to compute fast approximate exponentials
     */
    MultiShapeletMatrixBuilder(
        MultiShapeletBasis const & basis,
        MultiShapeletFunction const & psf,
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        bool useApproximateExp = false
    );

    /**
     *  @brief Evaluate the matrix given the ellipse parameters of the model
     *
     *  @param[out]  output   Matrix to fill, with dimensions (x.getSize<0>(), basis.getSize()).
     *                        Will be zeroed before filling.
     *  @param[in]   ellipse  Ellipse parameters of the model, with center relative to the x and y
     *                        arrays passed at construction.
     */
    void build(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const;

    // dtor needs to be explicit in cc file, because Impl dtor is undefined in .h file
    ~MultiShapeletMatrixBuilder();

private:
    class Impl;
    boost::scoped_ptr<Impl> _impl;
};

}} // namespace lsst::shapelet

#endif // !LSST_SHAPELET_MultiShapeletBasis_h_INCLUDED
