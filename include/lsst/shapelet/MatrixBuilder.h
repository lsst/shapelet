// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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
#ifndef LSST_SHAPELET_MatrixBuilder_h_INCLUDED
#define LSST_SHAPELET_MatrixBuilder_h_INCLUDED

#include "ndarray.h"

#include "lsst/shapelet/constants.h"
#include "lsst/shapelet/ShapeletFunction.h"
#include "lsst/shapelet/MultiShapeletFunction.h"

namespace lsst { namespace shapelet {

class MultiShapeletBasis;

template <typename T>
class MatrixBuilderFactory;

template <typename T>
class MatrixBuilderWorkspace;

/**
 *  @brief Class that evaluates a (multi-)shapelet basis at predefined points
 *
 *  The output "matrix" has pixels along the rows, and basis elements along columns;
 *  this is the design matrix involved in a linear least squares fit for the basis
 *  coefficients.
 *
 *  A MatrixBuilder can be constructed in two ways: via one of its own constructors,
 *  or via a MatrixBuilderFactory.  Using the latter allows the workspace arrays used
 *  by the MatrixBuilder to be shared between instances.  See MatrixBuilderFactory
 *  and MatrixBuilderWorkspace for more information.
 */
template <typename T>
class MatrixBuilder {
public:

    typedef MatrixBuilderFactory<T> Factory;     ///< Factory type associated with this builder
    typedef MatrixBuilderWorkspace<T> Workspace; ///< Workspace type associated with this builder

    class Impl; // Private implementation, defined in .cc; public to allow us to subclass it there

    /**
     *  Create a MatrixBuilder that evaluates a simple non-compound shapelet basis.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] order      order of the shapelet basis
     */
    MatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order
    );

    /**
     *  Create a MatrixBuilder that evaluates a simple non-compound shapelet basis after convolving it
     *  with a ShapeletFunction.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] order      order of the shapelet basis
     *  @param[in] psf        function to convolve the basis with
     */
    MatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order,
        ShapeletFunction const & psf
    );

    /**
     *  Create a MatrixBuilder that evaluates a MultiShapeletBasis object.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] basis      basis object defining the functions the matrix evaluates
     */
    MatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletBasis const & basis
    );

    /**
     *  Create a MatrixBuilder that evaluates a MultiShapeletBasis object after convolving it with
     *  a MultiShapeletFunction.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] basis      basis object defining the functions the matrix evaluates
     *  @param[in] psf        function to convolve the basis with
     */
    MatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletBasis const & basis,
        MultiShapeletFunction const & psf
    );

    /// Return the number of data points
    int getDataSize() const;

    /// Return the number of basis elements
    int getBasisSize() const;

    /// Return a matrix appropriate for use as an output for operator().
    ndarray::Array<T,2,-2> allocateOutput() const;

    /**
     *  @brief Fill an array with the model matrix.
     *
     *  @param[out]  output   Matrix to fill, with dimensions (getDataSize(), getBasisSize()).
     *                        Will be zeroed before filling.
     *  @param[in]   ellipse  Ellipse parameters of the model, with center relative to the x and y
     *                        arrays passed at construction.
     */
    void operator()(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const;

    /**
     *  @brief Return a newly-allocated model matrix.
     *
     *  @param[in]   ellipse  Ellipse parameters of the model, with center relative to the x and y
     *                        arrays passed at construction.
     */
    ndarray::Array<T,2,-2> operator()(afw::geom::ellipses::Ellipse const & ellipse) const {
        ndarray::Array<T,2,-2> output = allocateOutput();
        (*this)(output, ellipse);
        return output;
    }

private:

    template <typename U> friend class MatrixBuilderFactory;

    explicit MatrixBuilder(std::shared_ptr<Impl> impl);

    std::shared_ptr<Impl> _impl;
};

/**
 *  @brief Reusable, shareable workspace for MatrixBuilder.
 *
 *  Multiple MatrixBuilders are often used together to evaluate a multi-shapelet model,
 *  and in this case it's more efficient for the MatrixBuilders to use the same memory
 *  for temporary arrays rather than have them each allocate their own workspace.  At other
 *  times, it may be useful to use one workspace for a sequence of throwaway
 *  MatrixBuilders, to avoid unnecessary memory allocations.  This class manages the
 *  memory used for a MatrixBuilder's workspace arrays, and provides methods for tracking
 *  it and sharing it between multple MatrixBuilders.  See MatrixBuilderFactory for examples.
 *
 *  MatrixBuilderWorkspace holds a ndarray::Manager::Ptr that "owns" the block of memory.  This
 *  is passed on to any MatrixBuilders constructed with it, ensuring that the block of memory
 *  remains alive with the MatrixBuilder even if the workspace object goes out of scope - so it
 *  is not necessary to keep the workspace object alive manually for the duration of its users.
 *
 *  In addition, MatrixBuilderWorkspace tracks both the current point in memory, and increments
 *  this as workspace matrices and vectors are created from it.  It also checks that any created
 *  arrays do not exceed the bounds of the allocated memory.  In order to share workspace memory
 *  between MatrixBuilders, the workspace object is simply copied - copied objects share the
 *  same memory, but maintain different "current" pointers, and hence create arrays from the
 *  same block of memory.
 */
template <typename T>
class MatrixBuilderWorkspace {
public:

    typedef Eigen::Map< Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> > Matrix; ///< Workspace matrix type
    typedef Eigen::Map< Eigen::Array<T,Eigen::Dynamic,1> > Vector;  ///< Workspace vector type

    /// Allocate a block of memory with the given number of elements (of type T).
    explicit MatrixBuilderWorkspace(int size);

    /**
     *  Copy the current state of the workspace, allowing multiple MatrixBuilders to
     *  start their workspace arrays at the same point in memory and hence share workspace.
     *
     *  See MatrixBuilderFactory for an example.
     */
    MatrixBuilderWorkspace(MatrixBuilderWorkspace const & other) :
        _current(other._current),
        _end(other._end),
        _manager(other._manager)
    {}

    /// Return the size (in units of sizeof(T)) of the unused memory in the workspace.
    int getRemaining() const { return _end - _current; }

#ifndef SWIG

    /// Create a matrix from the workspace memory, and increment the current location accordingly.
    Matrix makeMatrix(int rows, int cols);

    /// Create a vector from the workspace memory, and increment the current location accordingly.
    Vector makeVector(int size);

    /// Manually increment the current location.
    void increment(int size);

    /// Return the manager object that owns the block of memory.
    ndarray::Manager::Ptr getManager() const { return _manager; }

#endif

private:

    void operator=(MatrixBuilderWorkspace & other); // disabled

    T * _current;
    T * _end;
    ndarray::Manager::Ptr _manager;
};

/**
 *  A factory class for MatrixBuilder, providing more control over workspace memory.
 *
 *  To allocate workspace manually for a MatrixBuilder, we use the following pattern:;
 *  @code
 *  MatrixBuilderFactory<T> factory(...);
 *  MatrixBuilderWorkspace<T> workspace(factory.computeWorkspace());
 *  MatrixBuilder<T> builder = factory(workspace);
 *  @endcode
 *  Because we aren't doing anything special with the workspace, however, this is actually
 *  exactly equivalent to just constructing the MatrixBuilder directly, i.e.:
 *  @code
 *  MatrixBuilder<T> builder(...);
 *  @endcode
 *
 *  A more interesting case is if we want to use the same workspace for a pair of MatrixBuilders:
 *  @code
 *  MatrixBuilderFactory<T> factory1(...);
 *  MatrixBuilderFactory<T> factory2(...);
 *  MatrixBuilderWorkspace<T> workspace1(
 *     std::max(factory1.computeWorkspace(), factory2.computeWorkspace())
 *  );
 *  MatrixBuilderWorkspace<T> workspace2(workspace1);
 *  MatrixBuilder<T> builder1 = factory1(workspace1);
 *  MatrixBuilder<T> builder2 = factory2(workspace2);
 *  @endcode
 */
template <typename T>
class MatrixBuilderFactory {
public:

    typedef MatrixBuilderWorkspace<T> Workspace; ///< Associated workspace class
    typedef MatrixBuilder<T> Builder; ///< Associated builder class

    class Impl; // Private implementation, defined in .cc; public to allow us to subclass it there

    /**
     *  Create a MatrixBuilder that evaluates a simple non-compound shapelet basis.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] order      order of the shapelet basis
     */
    MatrixBuilderFactory(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order
    );

    /**
     *  Create a MatrixBuilder that evaluates a simple non-compound shapelet basis after convolving it
     *  with a ShapeletFunction.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] order      order of the shapelet basis
     *  @param[in] psf        function to convolve the basis with
     */
    MatrixBuilderFactory(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order,
        ShapeletFunction const & psf
    );

    /**
     *  Create a MatrixBuilder that evaluates a MultiShapeletBasis object.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] basis      basis object defining the functions the matrix evaluates
     */
    MatrixBuilderFactory(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletBasis const & basis
    );

    /**
     *  Create a MatrixBuilder that evaluates a MultiShapeletBasis object after convolving it with
     *  a MultiShapeletFunction.
     *
     *  @param[in] x          column positions at which the basis should be evaluated.
     *  @param[in] y          row positions at which the basis should be evaluated (same size as x).
     *  @param[in] basis      basis object defining the functions the matrix evaluates
     *  @param[in] psf        function to convolve the basis with
     */
    MatrixBuilderFactory(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletBasis const & basis,
        MultiShapeletFunction const & psf
    );

    /// Return the number of data points
    int getDataSize() const;

    /// Return the number of basis elements
    int getBasisSize() const;

    /// Return the size of the workspace needed for this MatrixBuilder, in elements of T
    int computeWorkspace() const;

    /// Return a new MatrixBuilder with internal, unshared workspace
    MatrixBuilder<T> operator()() const;

    /// Return a new MatrixBuilder using the given workspace
    MatrixBuilder<T> operator()(Workspace & workspace) const;

private:
    std::shared_ptr<Impl> _impl;
};

}} // namespace lsst::shapelet

#endif // !LSST_SHAPELET_MatrixBuilder_h_INCLUDED
