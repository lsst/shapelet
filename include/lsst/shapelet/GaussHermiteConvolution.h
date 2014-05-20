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

#ifndef LSST_AFW_MATH_SHAPELETS_GaussHermiteConvolution
#define LSST_AFW_MATH_SHAPELETS_GaussHermiteConvolution

#include "ndarray.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/shapelet/HermiteTransformMatrix.h"

#include <boost/scoped_ptr.hpp>

namespace lsst { namespace shapelet {

class ShapeletFunction;

/**
 *  @brief A parametrized matrix that performs a convolution in shapelet space.
 *
 *  GaussHermiteConvolution is defined only for the HERMITE basis type.
 */
class GaussHermiteConvolution : private boost::noncopyable {
public:

    /**
     *  @brief Evaluate a shapelet convolution matrix in the given array.
     *
     *  @param[in,out] ellipse   On input, the ellipse core of the unconvolved shapelet expansion.
     *                           On output, the ellipse core of the convolved shapelet expansion.
     *
     *  The returned array is owned by the GaussHermiteConvolution object and will be modified
     *  the next time evaluate() is called.
     */
    ndarray::Array<double const,2,2> evaluate(afw::geom::ellipses::Ellipse & ellipse) const;

    /// @brief Return the order of the to-be-convolved shapelet basis.
    int getColOrder() const;

    /// @brief Return the order of the post-convolution shapelet basis.
    int getRowOrder() const;

    /// @brief Construct a matrix that convolves a basis of the given order with the given shapelet function.
    GaussHermiteConvolution(int colOrder, ShapeletFunction const & psf);

    // Must be defined in .cc file so it can see Impl dtor.
    ~GaussHermiteConvolution();

#ifndef SWIG
    class Impl;
#endif

private:

    boost::scoped_ptr<Impl> _impl;
};


}} // namespace lsst::shapelet

#endif // !LSST_AFW_MATH_SHAPELETS_GaussHermiteConvolution
