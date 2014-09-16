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

//===========================================================================================================
//================== Implementation Base Classes ============================================================
//===========================================================================================================

template <typename T>
class MatrixBuilder<T>::Impl {
public:

    virtual int getDataSize() const = 0;

    virtual int getBasisSize() const = 0;

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) = 0;

    virtual ~Impl() {}

};

template <typename T>
class MatrixBuilderFactory<T>::Impl {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;
    typedef typename MatrixBuilder<T>::Impl BuilderImpl;

    virtual int getDataSize() const = 0;

    virtual int getBasisSize() const = 0;

    virtual int computeWorkspace() const = 0;

    PTR(BuilderImpl) apply() const {
        Workspace workspace(computeWorkspace());
        return apply(workspace, workspace.getManager());
    }

    virtual PTR(BuilderImpl) apply(
        Workspace & workspace,
        ndarray::Manager::Ptr manager=ndarray::Manager::Ptr()
    ) const = 0;

    virtual ~Impl() {}

};

//===========================================================================================================
//================== Single-Component Implementation Base Classes ===========================================
//===========================================================================================================

namespace {

template <typename T>
class SimpleMatrixBuilderImpl : public MatrixBuilder<T>::Impl {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    SimpleMatrixBuilderImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        Workspace * workspace,
        ndarray::Manager::Ptr const & manager
    ) : _x(x), _y(y),
        _xt(workspace->makeVector(_x.template getSize<0>())),
        _yt(workspace->makeVector(_y.template getSize<0>())),
        _detFactor(1.0),
        _manager(manager)
    {}

    virtual int getDataSize() const { return _x.template getSize<0>(); }

    void readEllipse(afw::geom::ellipses::Ellipse const & ellipse) {
        afw::geom::AffineTransform transform = ellipse.getGridTransform();
        _xt = _x.template asEigen<Eigen::ArrayXpr>() * transform[afw::geom::AffineTransform::XX]
            + _y.template asEigen<Eigen::ArrayXpr>() * transform[afw::geom::AffineTransform::XY]
            + transform[afw::geom::AffineTransform::X];
        _yt = _x.template asEigen<Eigen::ArrayXpr>() * transform[afw::geom::AffineTransform::YX]
            + _y.template asEigen<Eigen::ArrayXpr>() * transform[afw::geom::AffineTransform::YY]
            + transform[afw::geom::AffineTransform::Y];
        _detFactor = transform.getLinear().computeDeterminant();
    }

private:
    ndarray::Array<T const,1,1> _x;
    ndarray::Array<T const,1,1> _y;
protected:
    Eigen::Map< Eigen::Array<T,Eigen::Dynamic,1> > _xt;
    Eigen::Map< Eigen::Array<T,Eigen::Dynamic,1> > _yt;
    T _detFactor;
private:
    ndarray::Manager::Ptr _manager;
};

template <typename T>
class SimpleMatrixBuilderFactoryImpl : public MatrixBuilderFactory<T>::Impl {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    virtual int getDataSize() const { return _x.template getSize<0>(); }

    virtual int computeWorkspace() const { return 2*_x.template getSize<0>(); }

    SimpleMatrixBuilderFactoryImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y
    ) : _x(x), _y(y) {}

protected:
    ndarray::Array<T const,1,1> _x;
    ndarray::Array<T const,1,1> _y;
};

} // anonymous

//===========================================================================================================
//================== Non-Convolved, Non-Remapped Shapelet Implementation ====================================
//===========================================================================================================

namespace {

template <typename T>
class ShapeletMatrixBuilderImpl : public SimpleMatrixBuilderImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    ShapeletMatrixBuilderImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order,
        Workspace * workspace,
        ndarray::Manager::Ptr const & manager
    ) : SimpleMatrixBuilderImpl<T>(x, y, workspace, manager),
        _order(order),
        _gaussian(workspace->makeVector(x.template getSize<0>())),
        _xHermite(workspace->makeMatrix(x.template getSize<0>(), order + 1)),
        _yHermite(workspace->makeMatrix(x.template getSize<0>(), order + 1))
    {}

    virtual int getBasisSize() const { return computeSize(_order); }

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        apply(output.template asEigen<Eigen::ArrayXpr>(), ellipse);
    }

    template <typename EigenArrayT>
    void apply(
        EigenArrayT output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        this->readEllipse(ellipse);
        fillGaussian();
        fillHermite1d(_xHermite, this->_xt);
        fillHermite1d(_yHermite, this->_yt);
        for (PackedIndex i; i.getOrder() <= _order; ++i) {
            output.col(i.getIndex())
                += this->_detFactor * _gaussian * _xHermite.col(i.getX()) * _yHermite.col(i.getY());
        }
    }

    void fillGaussian() {
        _gaussian = (-0.5*(this->_xt.square() + this->_yt.square())).exp();
    }

    template <typename CoordArray>
    void fillHermite1d(
        typename Workspace::Matrix & output,
        CoordArray const & coord
    ) {
        if (output.cols() > 0) {
            output.col(0).setConstant(BASIS_NORMALIZATION);
        }
        if (output.cols() > 1) {
            output.col(1) = intSqrt(2) * coord * output.col(0);
        }
        for (int j = 2; j <= _order; ++j) {
            output.col(j) = rationalSqrt(2, j) * coord * output.col(j-1)
                - rationalSqrt(j - 1, j) * output.col(j-2);
        }
    }

protected:
    int _order;
    typename Workspace::Vector _gaussian;
    typename Workspace::Matrix _xHermite;
    typename Workspace::Matrix _yHermite;
};

template <typename T>
class ShapeletMatrixBuilderFactoryImpl : public SimpleMatrixBuilderFactoryImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;
    typedef ShapeletMatrixBuilderImpl<T> BuilderImpl;

    ShapeletMatrixBuilderFactoryImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order
    ) : SimpleMatrixBuilderFactoryImpl<T>(x, y),
        _order(order),
        _basisSize(computeSize(order))
    {}

    virtual int getBasisSize() const { return _basisSize; }

    virtual int computeWorkspace() const {
        return this->getDataSize()*(1 + 2*(_order + 1))
            + SimpleMatrixBuilderFactoryImpl<T>::computeWorkspace();
    }

    virtual PTR(typename MatrixBuilder<T>::Impl) apply(
        Workspace & workspace,
        ndarray::Manager::Ptr manager=ndarray::Manager::Ptr()
    ) const {
        return boost::make_shared<BuilderImpl>(this->_x, this->_y, _order, &workspace, manager);
    }

protected:
    int _order;
    int _basisSize;
};

} // anonymous


//===========================================================================================================
//================== Non-Convolved, Remapped Shapelet Implementation ========================================
//===========================================================================================================

namespace {

template <typename T>
class RemappedShapeletMatrixBuilderImpl : public ShapeletMatrixBuilderImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    RemappedShapeletMatrixBuilderImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order,
        double radius,
        ndarray::Array<double const,2,2> const & remapMatrix,
        Workspace * workspace,
        ndarray::Manager::Ptr const & manager
    ) : ShapeletMatrixBuilderImpl<T>(x, y, order, workspace, manager),
        _radius(radius),
        _ellipse(afw::geom::ellipses::Quadrupole()),
        _remapMatrix(remapMatrix.asEigen().cast<T>().transpose()), // transpose to preserve memory order
        _lhs(workspace->makeMatrix(x.template getSize<0>(), computeSize(order)))
    {}

    virtual int getBasisSize() const { return _remapMatrix.rows(); }

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        _ellipse = ellipse;
        _ellipse.scale(_radius);
        _lhs.setZero();
        ShapeletMatrixBuilderImpl<T>::apply(_lhs, _ellipse);
        output.asEigen() += _lhs.matrix() * _remapMatrix.transpose(); // undo the transpose in the ctor
    }

protected:
    double _radius;
    mutable afw::geom::ellipses::Ellipse _ellipse;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> _remapMatrix;
    typename Workspace::Matrix _lhs;
};

template <typename T>
class RemappedShapeletMatrixBuilderFactoryImpl : public ShapeletMatrixBuilderFactoryImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;
    typedef RemappedShapeletMatrixBuilderImpl<T> BuilderImpl;

    RemappedShapeletMatrixBuilderFactoryImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        int order,
        double radius,
        ndarray::Array<double const,2,2> const & remapMatrix
    ) : ShapeletMatrixBuilderFactoryImpl<T>(x, y, order),
        _radius(radius),
        _remapMatrix(remapMatrix)
    {}

    virtual int getBasisSize() const { return _remapMatrix.getSize<1>(); }

    virtual int computeWorkspace() const {
        return ShapeletMatrixBuilderFactoryImpl<T>::computeWorkspace()
            + ShapeletMatrixBuilderFactoryImpl<T>::getBasisSize() * this->getDataSize();
    }

    virtual PTR(typename MatrixBuilder<T>::Impl) apply(
        Workspace & workspace,
        ndarray::Manager::Ptr manager=ndarray::Manager::Ptr()
    ) const {
        return boost::make_shared<BuilderImpl>(
            this->_x, this->_y, this->_order, _radius, _remapMatrix, &workspace, manager
        );
    }

protected:
    double _radius;
    ndarray::Array<double const,2,2> _remapMatrix;
};

} // anonymous


//===========================================================================================================
//================== Multi-Component Implementations ========================================================
//===========================================================================================================

namespace {

template <typename T>
class CompoundMatrixBuilderImpl : public MatrixBuilder<T>::Impl {
public:

    typedef typename MatrixBuilder<T>::Impl Component;
    typedef std::vector<PTR(Component)> Vector;
    typedef typename Vector::const_iterator Iterator;

    explicit CompoundMatrixBuilderImpl(Vector const & components) : _components(components) {}

    virtual int getDataSize() const { return _components.front()->getDataSize(); }

    virtual int getBasisSize() const { return _components.front()->getBasisSize(); }

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        for (Iterator i = _components.begin(); i != _components.end(); ++i) {
            (**i).apply(output, ellipse);
        }
    }

private:
    Vector _components;
};

template <typename T>
class CompoundMatrixBuilderFactoryImpl : public MatrixBuilderFactory<T>::Impl {
public:

    typedef typename MatrixBuilderFactory<T>::Impl Component;
    typedef std::vector<PTR(Component)> Vector;
    typedef typename Vector::const_iterator Iterator;

    typedef MatrixBuilderWorkspace<T> Workspace;
    typedef CompoundMatrixBuilderImpl<T> BuilderImpl;

    CompoundMatrixBuilderFactoryImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletBasis const & basis
    ) {
        _components.reserve(basis.getComponentCount());
        for (MultiShapeletBasis::Iterator i = basis.begin(); i != basis.end(); ++i) {
            _components.push_back(
                boost::make_shared< RemappedShapeletMatrixBuilderFactoryImpl<T> >(
                    x, y, i->getOrder(), i->getRadius(), i->getMatrix()
                )
            );
        }
    }

    CompoundMatrixBuilderFactoryImpl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        MultiShapeletBasis const & basis,
        MultiShapeletFunction const & psf
    ) {
        _components.reserve(psf.getElements().size() * basis.getComponentCount());
        // TODO
    }

    virtual int getBasisSize() const { return _components.front()->getBasisSize(); }

    virtual int getDataSize() const { return _components.front()->getDataSize(); }

    virtual int computeWorkspace() const {
        int ws = 0;
        for (Iterator i = _components.begin(); i != _components.end(); ++i) {
            ws = std::max((**i).computeWorkspace(), ws);
        }
        return ws;
    }

    virtual PTR(typename MatrixBuilder<T>::Impl) apply(
        Workspace & workspace,
        ndarray::Manager::Ptr manager=ndarray::Manager::Ptr()
    ) const {
        typename BuilderImpl::Vector builders;
        builders.reserve(_components.size());
        for (Iterator i = _components.begin(); i != _components.end(); ++i) {
            // By copying the workspace here, we prevent the original from having its pointer updated (yet),
            // and hence each component builder starts grabbing workspace arrays from the same point,
            // and they all end up sharing the same space.  That's what we want, because we call them one
            // at a time, and they don't need anything to remain between calls.
            Workspace wsCopy(workspace);
            builders.push_back((**i).apply(wsCopy, manager));
        }
        // Now, at the end, we increment the workspace by the amount we claimed we needed, which was
        // the maximum needed by any individual components.
        workspace.increment(computeWorkspace());
        return boost::make_shared<BuilderImpl>(builders);
    }

private:
    Vector _components;
};

} // anonymous

//===========================================================================================================
//================== MatrixBuilder ==========================================================================
//===========================================================================================================

template <typename T>
int MatrixBuilder<T>::getDataSize() const {
    return _impl->getDataSize();
}

template <typename T>
int MatrixBuilder<T>::getBasisSize() const {
    return _impl->getBasisSize();
}

template <typename T>
ndarray::Array<T,2,-2> MatrixBuilder<T>::allocateOutput() const {
    ndarray::Array<T,2,2> t = ndarray::allocate(getBasisSize(), getDataSize());
    t.deep() = 0.0;
    return t.transpose();
}

template <typename T>
void MatrixBuilder<T>::operator()(
    ndarray::Array<T,2,-1> const & output,
    afw::geom::ellipses::Ellipse const & ellipse
) const {
    _impl->apply(output, ellipse);
}

template <typename T>
MatrixBuilder<T>::MatrixBuilder(PTR(Impl) impl) :
    _impl(impl)
{}


//===========================================================================================================
//================== MatrixBuilderWorkspace =================================================================
//===========================================================================================================

template <typename T>
MatrixBuilderWorkspace<T>::MatrixBuilderWorkspace(int size) {
    std::pair<ndarray::Manager::Ptr,T*> pair = ndarray::SimpleManager<T>::allocate(size);
    _current = pair.second;
    _end = pair.second + size;
    _manager = pair.first;
}


template <typename T>
typename MatrixBuilderWorkspace<T>::Matrix MatrixBuilderWorkspace<T>::makeMatrix(int rows, int cols) {
    Matrix m(_current, rows, cols);
    _current += rows*cols;
    if (_current > _end) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            "Allocated workspace is too small"
        );
    }
    return m;
}

template <typename T>
typename MatrixBuilderWorkspace<T>::Vector MatrixBuilderWorkspace<T>::makeVector(int size) {
    Vector v(_current, size);
    _current += size;
    if (_current > _end) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            "Allocated workspace is too small"
        );
    }
    return v;
}

template <typename T>
void MatrixBuilderWorkspace<T>::increment(int size) {
    _current += size;
}

//===========================================================================================================
//================== MatrixBuilderFactory ===================================================================
//===========================================================================================================

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order
) : _impl(boost::make_shared< ShapeletMatrixBuilderFactoryImpl<T> >(x, y, order))
{}

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order,
    ShapeletFunction const & psf
) {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not implemented");
}

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order,
    MultiShapeletFunction const & psf
) {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not implemented");
}

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletBasis const & basis
) {
    if (basis.getComponentCount() == 1) {
        MultiShapeletBasisComponent const & component = *basis.begin();
        ndarray::EigenView<double const,2,2> matrix(component.getMatrix());
        if (
            std::abs(component.getRadius() - 1.0) < std::numeric_limits<double>::epsilon() &&
            matrix.rows() == matrix.cols() &&
            matrix.isApprox(
                Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()),
                std::numeric_limits<double>::epsilon()
            )
        ) {
            _impl = boost::make_shared< ShapeletMatrixBuilderFactoryImpl<T> >(x, y, component.getOrder());
        } else {
            _impl = boost::make_shared< RemappedShapeletMatrixBuilderFactoryImpl<T> >(
                x, y, component.getOrder(), component.getRadius(), component.getMatrix()
            );
        }
    } else {
        _impl = boost::make_shared< CompoundMatrixBuilderFactoryImpl<T> >(x, y, basis);
    }
}

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletFunction const & psf,
    MultiShapeletBasis const & basis
) {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not implemented");
}

template <typename T>
int MatrixBuilderFactory<T>::getDataSize() const { return _impl->getDataSize(); }

template <typename T>
int MatrixBuilderFactory<T>::getBasisSize() const { return _impl->getBasisSize(); }

template <typename T>
int MatrixBuilderFactory<T>::computeWorkspace() const { return _impl->computeWorkspace(); }

template <typename T>
MatrixBuilder<T> MatrixBuilderFactory<T>::operator()() const {
    return MatrixBuilder<T>(_impl->apply());
}

template <typename T>
MatrixBuilder<T> MatrixBuilderFactory<T>::operator()(Workspace & workspace) const {
    return MatrixBuilder<T>(_impl->apply(workspace));
}

//===========================================================================================================
//================== Explicit Instantiation =================================================================
//===========================================================================================================

#define INSTANTIATE(T)                                          \
    template class MatrixBuilder<T>;                            \
    template class MatrixBuilderFactory<T>;                     \
    template class MatrixBuilderWorkspace<T>;                   \
    template class SimpleMatrixBuilderImpl<T>;                  \
    template class SimpleMatrixBuilderFactoryImpl<T>;           \
    template class ShapeletMatrixBuilderImpl<T>;                \
    template class ShapeletMatrixBuilderFactoryImpl<T>;         \
    template class RemappedShapeletMatrixBuilderImpl<T>;        \
    template class RemappedShapeletMatrixBuilderFactoryImpl<T>; \
    template class CompoundMatrixBuilderImpl<T>;                \
    template class CompoundMatrixBuilderFactoryImpl<T>

INSTANTIATE(float);
INSTANTIATE(double);

}} // namespace lsst::shapelet
