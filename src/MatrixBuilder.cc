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

#include <cmath>
#include "boost/make_shared.hpp"
#include "ndarray/eigen.h"

#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/shapelet/GaussHermiteConvolution.h"

namespace lsst { namespace shapelet {

namespace {

template <typename T>
class EllipseHelper {
public:

    explicit EllipseHelper(int dataSize) : detFactor(1.0), xyt(dataSize, 0) {}

    void readEllipse(
        Eigen::Matrix<T,Eigen::Dynamic,2> const & xy,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        afw::geom::AffineTransform transform = ellipse.getGridTransform();
        xyt.transpose() = xy.transpose() * transform.getLinear().getMatrix().transpose().template cast<T>();
        xyt.rowwise() += transform.getTranslation().asEigen().transpose().template cast<T>();
        detFactor = transform.getLinear().computeDeterminant();
    }

    T detFactor;
    Eigen::Matrix<T,Eigen::Dynamic,2> xyt;
};

template <typename T>
class GaussianHelper {
public:

    void apply(
        EllipseHelper<T> const & ellipseHelper,
        ndarray::Array<T,1,1> const & output
    ) {
        static T const NORM = 1.0 / std::sqrt(M_PI); // normalization to match shapelets
        // TODO: check that rowwise().squaredNorm() is optimized as well as explicitly writing it as
        // coeffwise array operations (may depend on whether we transpose xy).
        output.template asEigen<Eigen::ArrayXpr>() +=
            (-0.5*ellipseHelper.xyt.rowwise().squaredNorm()).array().exp()
            * ellipseHelper.detFactor
            * NORM;
    }

};

template <typename T>
class ShapeletHelper {
public:

    explicit ShapeletHelper(int order, int dataSize) :
        _order(order),
        _expWorkspace(dataSize),
        _xWorkspace(order + 1, dataSize),
        _yWorkspace(order + 1, dataSize)
    {}

    void apply(
        EllipseHelper<T> const & ellipseHelper,
        ndarray::Array<T,2,-1> const & output
    ) {
        _expWorkspace =
            (-0.5*ellipseHelper.xyt.rowwise().squaredNorm()).array().exp() * ellipseHelper.detFactor;
        _fillHermite1d(_xWorkspace, ellipseHelper.xyt.col(0).array());
        _fillHermite1d(_yWorkspace, ellipseHelper.xyt.col(1).array());
        ndarray::EigenView<T,2,-1,Eigen::ArrayXpr> view(output);
        for (PackedIndex i; i.getOrder() <= _order; ++i) {
            view.col(i.getIndex()) += _expWorkspace*_xWorkspace.col(i.getX())*_yWorkspace.col(i.getY());
        }
    }

private:

    template <typename CoordArray>
    void _fillHermite1d(
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & workspace,
        CoordArray const & coord
    ) {
        if (_order >= workspace.cols()) {
            workspace.resize(coord.size(), _order + 1);
        }
        if (workspace.cols() > 0) {
            workspace.col(0).setConstant(BASIS_NORMALIZATION);
        }
        if (workspace.cols() > 1) {
            workspace.col(1) = intSqrt(2) * coord * workspace.col(0);
        }
        for (int j = 2; j <= _order; ++j) {
            workspace.col(j) = rationalSqrt(2, j) * coord * workspace.col(j-1)
                - rationalSqrt(j - 1, j) * workspace.col(j-2);
        }
    }

    int _order;
    Eigen::Array<T,Eigen::Dynamic,1> _expWorkspace;
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> _xWorkspace;
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> _yWorkspace;
};

template <typename T>
class GaussianMatrixBuilder : public MatrixBuilder<T> {
public:

    GaussianMatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y
    ) : MatrixBuilder<T>(x, y, 1),
        _ellipseHelper(this->getDataSize())
    {}

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        output.deep() = 0.0;
        _ellipseHelper.readEllipse(this->_xy, ellipse);
        _gaussianHelper.apply(_ellipseHelper, output.transpose()[0]);
    }

private:
    mutable EllipseHelper<T> _ellipseHelper;
    mutable GaussianHelper<T> _gaussianHelper;
};

template <typename T>
class ConvolvedGaussianMatrixBuilder : public MatrixBuilder<T> {
public:

    ConvolvedGaussianMatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        afw::geom::ellipses::Ellipse const & psfEllipse,
        double psfCoefficient
    ) : MatrixBuilder<T>(x, y, 1),
        _ellipseHelper(this->getDataSize()),
        _psfEllipse(psfEllipse),
        _psfCoefficient(psfCoefficient)
    {}

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        output.deep() = 0.0;
        _ellipseHelper.readEllipse(this->_xy, ellipse.convolve(_psfEllipse));
        _gaussianHelper.apply(_ellipseHelper, output.transpose()[0]);
        output.asEigen() *= _psfCoefficient / ShapeletFunction::FLUX_FACTOR;
    }

private:
    mutable EllipseHelper<T> _ellipseHelper;
    mutable GaussianHelper<T> _gaussianHelper;
    afw::geom::ellipses::Ellipse _psfEllipse;
    double _psfCoefficient;
};

template <typename T>
class ShapeletMatrixBuilder : public MatrixBuilder<T> {
public:

    ShapeletMatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order
    ) : MatrixBuilder<T>(x, y, computeSize(order)),
        _ellipseHelper(this->getDataSize()),
        _shapeletHelper(order, this->getDataSize())
    {}

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        output.deep() = 0.0;
        _ellipseHelper.readEllipse(this->_xy, ellipse);
        _shapeletHelper.apply(_ellipseHelper, output);
    }

private:
    mutable EllipseHelper<T> _ellipseHelper;
    mutable ShapeletHelper<T> _shapeletHelper;
};


template <typename T>
class ConvolvedShapeletMatrixBuilder : public MatrixBuilder<T> {
public:

    ConvolvedShapeletMatrixBuilder(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        ShapeletFunction const & psf,
        int order
    ) : MatrixBuilder<T>(x, y, computeSize(order)),
        _convolution(order, psf),
        _convolutionWorkspace(this->getDataSize(), _convolution.getRowOrder()),
        _ellipseHelper(this->getDataSize()),
        _shapeletHelper(_convolution.getRowOrder(), this->getDataSize())
    {}

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        _convolutionWorkspace.deep() = 0.0;
        afw::geom::ellipses::Ellipse convolvedEllipse(ellipse);
        ndarray::Array<double const,2,2> convolutionMatrix = _convolution.evaluate(convolvedEllipse);
        _ellipseHelper.readEllipse(this->_xy, convolvedEllipse);
        _shapeletHelper.apply(_ellipseHelper, _convolutionWorkspace);
        output.asEigen() = _convolutionWorkspace.asEigen() * convolutionMatrix.asEigen().cast<T>();
    }

private:
    GaussHermiteConvolution _convolution;
    ndarray::Array<T,2,-1> _convolutionWorkspace;
    mutable EllipseHelper<T> _ellipseHelper;
    mutable ShapeletHelper<T> _shapeletHelper;
};

} // anonymous

template <typename T>
PTR(MatrixBuilder<T>) makeMatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order
) {
    if (order == 0) {
        return boost::make_shared< GaussianMatrixBuilder<T> >(x, y);
    } else {
        return boost::make_shared< ShapeletMatrixBuilder<T> >(x, y, order);
    }
}
template <typename T>
PTR(MatrixBuilder<T>) makeMatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    ShapeletFunction const & psf,
    int order
) {
    if (order == 0 && psf.getOrder() == 0) {
        return boost::make_shared< ConvolvedGaussianMatrixBuilder<T> >(
            x, y, psf.getEllipse(), psf.getCoefficients()[0]
        );
    } else {
        return boost::make_shared< ConvolvedShapeletMatrixBuilder<T> >(x, y, psf, order);
    }
}

template <typename T>
PTR(MatrixBuilder<T>) makeMatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletFunction const & psf,
    int order
) {
    if (psf.getElements().size() == 1u) {
        return makeMatrixBuilder(x, y, psf.getElements().front(), order);
    } else {
        throw pex::exceptions::LogicError("Not implemented");
    }
}

template <typename T>
PTR(MatrixBuilder<T>) makeMatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletBasis const & basis
) {
    throw pex::exceptions::LogicError("Not implemented");
}

template <typename T>
PTR(MatrixBuilder<T>) makeMatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletFunction const & psf,
    MultiShapeletBasis const & basis
) {
    throw pex::exceptions::LogicError("Not implemented");
}

#define INSTANTIATE(T)                                      \
    template class MatrixBuilder<T>;                        \
    template PTR(MatrixBuilder<T>) makeMatrixBuilder(       \
        ndarray::Array<T const,1,1> const &,                \
        ndarray::Array<T const,1,1> const &,                \
        int order                                           \
    );                                                      \
    PTR(MatrixBuilder<T>) makeMatrixBuilder(                \
        ndarray::Array<T const,1,1> const &,                \
        ndarray::Array<T const,1,1> const &,                \
        ShapeletFunction const &,                           \
        int order                                           \
    );                                                      \
    PTR(MatrixBuilder<T>) makeMatrixBuilder(                \
        ndarray::Array<T const,1,1> const &,                \
        ndarray::Array<T const,1,1> const &,                \
        MultiShapeletFunction const &,                      \
        int order                                           \
    );                                                      \
    PTR(MatrixBuilder<T>) makeMatrixBuilder(                \
        ndarray::Array<T const,1,1> const &,                \
        ndarray::Array<T const,1,1> const &,                \
        MultiShapeletBasis const &                          \
    );                                                      \
    PTR(MatrixBuilder<T>) makeMatrixBuilder(                \
        ndarray::Array<T const,1,1> const &,                \
        ndarray::Array<T const,1,1> const &,                \
        MultiShapeletFunction const &,                      \
        MultiShapeletBasis const &                          \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}} // namespace lsst::shapelet
