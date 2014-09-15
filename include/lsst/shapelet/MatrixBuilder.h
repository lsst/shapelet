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

/**
 *  @brief Class that evaluates a (multi-)shapelet basis at predefined points
 *
 *  The output "matrix" has pixels along the rows, and basis elements along columns;
 *  this is the design matrix involved in a linear least squares fit for the basis
 *  coefficients.
 */
template <typename T>
class MatrixBuilder {
public:

    class Impl;

    /// Return the number of data points
    int getDataSize() const;

    /// Return the number of basis elements
    int getBasisSize() const;

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

private:

    template <typename U> friend class MatrixBuilderFactory;

    explicit MatrixBuilder(PTR(Impl) impl);

    PTR(Impl) _impl;
};

template <typename T>
class MatrixBuilderWorkspace {
public:

    typedef Eigen::Map< Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> > Matrix;
    typedef Eigen::Map< Eigen::Array<T,Eigen::Dynamic,1> > Vector;

    explicit MatrixBuilderWorkspace(int size);

    bool check(T const * end) const { return end == _end; }

    Matrix makeMatrix(int rows, int cols) {
        Matrix m(_current, rows, cols);
        _current += rows*cols;
        return m;
    }

    Vector makeVector(int size) {
        Vector v(_current, size);
        _current += size;
        return v;
    }

    ndarray::Manager::Ptr getManager() const { return _manager; }

private:
    T * _current;
    T * _end;
    ndarray::Manager::Ptr _manager;
};

template <typename T>
class MatrixBuilderFactory {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    class Impl;

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
     *  Create a MatrixBuilder that evaluates a simple non-compound shapelet basis after convolving it
     *  with a MultiShapeletFunction.
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
        MultiShapeletFunction const & psf
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
     *  @param[in] psf        function to convolve the basis with
     *  @param[in] basis      basis object defining the functions the matrix evaluates
     */
    MatrixBuilderFactory(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletFunction const & psf,
        MultiShapeletBasis const & basis
    );

    /// Return the number of data points
    int getDataSize() const;

    /// Return the number of basis elements
    int getBasisSize() const;

    int computeWorkspace() const;

    MatrixBuilder<T> operator()() const;

    MatrixBuilder<T> operator()(Workspace & workspace) const;

private:
    PTR(Impl) _impl;
};

}} // namespace lsst::shapelet

#endif // !LSST_SHAPELET_MatrixBuilder_h_INCLUDED
