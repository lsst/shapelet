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
    ) const = 0;

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

    PTR(BuilderImpl) operator()() const {
        Workspace workspace(computeWorkspace());
        return (*this)(workspace, workspace.getManager());
    }

    virtual PTR(BuilderImpl) operator()(
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
        _expWorkspace(workspace->makeVector(x.template getSize<0>())),
        _xWorkspace(workspace->makeMatrix(x.template getSize<0>(), order + 1)),
        _yWorkspace(workspace->makeMatrix(y.template getSize<0>(), order + 1))
    {}

    virtual int getBasisSize() const { return computeSize(_order); }

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        // TODO
    }

private:
    int _order;
    typename Workspace::Vector _expWorkspace;
    typename Workspace::Matrix _xWorkspace;
    typename Workspace::Matrix _yWorkspace;
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
    ) : SimpleMatrixBuilderFactoryImpl<T>(x, y), _order(order), _basisSize(computeSize(order)) {}

    virtual int getBasisSize() const { return _basisSize; }

    virtual int computeWorkspace() const { return this->getDataSize()*(1 + 2*(_order + 1)); }

    virtual PTR(typename MatrixBuilder<T>::Impl) operator()(
        Workspace & workspace,
        ndarray::Manager::Ptr manager=ndarray::Manager::Ptr()
    ) const {
        return boost::make_shared<BuilderImpl>(this->_x, this->_y, _order, &workspace, manager);
    }

private:
    int _order;
    int _basisSize;
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

    CompoundMatrixBuilderImpl(int dataSize, int basisSize, Vector const & components) :
        _dataSize(dataSize), _basisSize(basisSize), _components(components)
    {}

    virtual int getDataSize() const { return _dataSize; }

    virtual int getBasisSize() const { return _basisSize; }

    virtual void apply(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        for (Iterator i = _components.begin(); i != _components.end(); ++i) {
            (**i).apply(output, ellipse);
        }
    }

private:
    int _dataSize;
    int _basisSize;
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


//===========================================================================================================
//================== Explicit Instantiation =================================================================
//===========================================================================================================

#define INSTANTIATE(T)                                  \
    template class MatrixBuilder<T>;                    \
    template class MatrixBuilderFactory<T>;             \
    template class SimpleMatrixBuilderImpl<T>;          \
    template class SimpleMatrixBuilderFactoryImpl<T>;   \
    template class ShapeletMatrixBuilderImpl<T>;        \
    template class ShapeletMatrixBuilderFactoryImpl<T>

INSTANTIATE(float);
INSTANTIATE(double);

}} // namespace lsst::shapelet
