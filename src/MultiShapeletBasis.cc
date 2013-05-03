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

#include "ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/shapelet/GaussHermiteConvolution.h"

namespace lsst { namespace shapelet {

MultiShapeletBasisComponent::MultiShapeletBasisComponent(
    double radius,
    int order,
    ndarray::Array<double const,2,2> const & matrix
) : _radius(radius), _order(order), _matrix(ndarray::copy(matrix)) {
    if (_matrix.getSize<0>() != computeSize(order)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("MultiShapeletBasisComponent matrix has %d rows; expected %d for order=%d")
             % _matrix.getSize<0>() % computeSize(order) % order).str()
        );
    }
}

MultiShapeletBasis::MultiShapeletBasis(int size) : _size(size), _components() {}

void MultiShapeletBasis::addComponent(
    double radius,
    int order,
    ndarray::Array<double const,2,2> const & matrix
) {
    if (matrix.getSize<1>() != _size) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Component matrix has %d columns; basis size is %d")
             % matrix.getSize<1>() % _size).str()
        );
    }
    _components.push_back(Component(radius, order, matrix));
}

void MultiShapeletBasis::normalize() {
    Eigen::ArrayXd totals = Eigen::VectorXd::Zero(_size);
    for (Iterator i = begin(); i != end(); ++i) {
        GaussHermiteEvaluator ev(i->getOrder());
        for (int n = 0; n < _size; ++n) {
            totals[n] += ev.sumIntegration(i->getMatrix()[ndarray::view()(n)]);
        }
    }
    for (ComponentVector::iterator i = _components.begin(); i != _components.end(); ++i) {
        ndarray::Array<double,2,2> newMatrix = ndarray::copy(i->getMatrix());
        for (int n = 0; n < _size; ++n) {
            newMatrix.asEigen().col(n) /= totals[n];
        }
        i->_matrix = newMatrix;
    }
}

void MultiShapeletBasis::scale(double factor) {
    for (ComponentVector::iterator i = _components.begin(); i != _components.end(); ++i) {
        i->_radius *= factor;
    }
}

void MultiShapeletBasis::merge(MultiShapeletBasis const & other) {
    ComponentVector newComponents;
    int newSize = _size + other.getSize();
    for (Iterator i = begin(); i != end(); ++i) {
        ndarray::Array<double,2,2> newMatrix = ndarray::allocate(i->getMatrix().getSize<0>(), newSize);
        newMatrix[ndarray::view()(0,_size)] = i->getMatrix();
        newMatrix[ndarray::view()(_size, _size + other.getSize())] = 0.0;
        newComponents.push_back(Component(i->getRadius(), i->getOrder(), newMatrix));
    }
    for (Iterator i = other.begin(); i != other.end(); ++i) {
        ndarray::Array<double,2,2> newMatrix = ndarray::allocate(i->getMatrix().getSize<0>(), newSize);
        newMatrix[ndarray::view()(0,_size)] = 0.0;
        newMatrix[ndarray::view()(_size, _size + other.getSize())] = i->getMatrix();
        newComponents.push_back(Component(i->getRadius(), i->getOrder(), newMatrix));
    }
    _size = newSize;
    _components.swap(newComponents);
}

MultiShapeletFunction MultiShapeletBasis::makeFunction(
    afw::geom::ellipses::Ellipse const & ellipse,
    ndarray::Array<double const,1,1> const & coefficients
) const {
    MultiShapeletFunction result;
    for (Iterator i = begin(); i != end(); ++i) {
        result.getElements().push_back(ShapeletFunction(i->getOrder(), HERMITE, ellipse));
        result.getElements().back().getEllipse().getCore().scale(i->getRadius());
        result.getElements().back().getCoefficients().asEigen()
            = i->getMatrix().asEigen() * coefficients.asEigen();
    }
    return result;
}

namespace {

template <typename T>
struct MatrixBuilderItem : public MultiShapeletBasis::Component {

    MatrixBuilderItem(
        MultiShapeletBasis::Component const & basisComponent,
        PTR(GaussHermiteConvolution) convolution_
    ) : MultiShapeletBasis::Component(basisComponent),
        convolution(convolution_)
    {}

    PTR(GaussHermiteConvolution) convolution;
    ndarray::Array<T,2,-1> workspace;
};

} // anonymous

template <typename T>
class MultiShapeletMatrixBuilder<T>::Impl {
public:

    typedef std::vector< MatrixBuilderItem<T> > ItemVector;

    Impl(
        ndarray::Array<T const,1,1> const & x,
        ndarray::Array<T const,1,1> const & y,
        bool useApproximateExp
    ) : modelBuilder(x, y, useApproximateExp) {}

    ModelBuilder<T> modelBuilder;
    ItemVector items;
};

template <typename T>
void MultiShapeletMatrixBuilder<T>::build(
    ndarray::Array<T,2,-1> const & output,
    afw::geom::ellipses::Ellipse const & ellipse
) const {
    output.asEigen().setZero();
    for (typename Impl::ItemVector::const_iterator i = _impl->items.begin(); i != _impl->items.end(); ++i) {
        afw::geom::ellipses::Ellipse itemEllipse(ellipse);
        itemEllipse.getCore().scale(i->getRadius());
        ndarray::Array<double const,2,2> convolution = i->convolution->evaluate(itemEllipse);
        i->workspace.asEigen().setZero();
        _impl->modelBuilder.update(itemEllipse);
        _impl->modelBuilder.addModelMatrix(i->convolution->getRowOrder(), i->workspace);
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rhs
            = (convolution.asEigen() * i->getMatrix().asEigen()).template cast<T>();
        output.asEigen() += i->workspace.asEigen() * rhs;
    }
}

template <typename T>
MultiShapeletMatrixBuilder<T>::MultiShapeletMatrixBuilder(
    MultiShapeletBasis const & basis,
    MultiShapeletFunction const & psf,
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    bool useApproximateExp
) : _impl(new Impl(x, y, useApproximateExp))
{
    // Add the cartesian product of (basis components) x (PSF elements) to the ItemVector
    int maxConvolvedOrder = 0;
    for (
        MultiShapeletFunction::ElementList::const_iterator psfIter = psf.getElements().begin();
        psfIter != psf.getElements().end();
        ++psfIter
    ) {
        int lastBasisOrder = -1;
        PTR(GaussHermiteConvolution) convolution;
        for (MultiShapeletBasis::Iterator basisIter = basis.begin(); basisIter != basis.end(); ++basisIter) {
            if (basisIter->getOrder() != lastBasisOrder) {
                convolution.reset(new GaussHermiteConvolution(basisIter->getOrder(), *psfIter));
                lastBasisOrder = basisIter->getOrder();
                maxConvolvedOrder = std::max(convolution->getRowOrder(), maxConvolvedOrder);
            }
            _impl->items.push_back(MatrixBuilderItem<T>(*basisIter, convolution));
        }
    }

    // Loop over the ItemVector one more time to set up workspace arrays with which to call
    // ModelBuilder.  We allocate one big array with the maximum workspace, then give each
    // item a view into it that only contains as much as it needs.
    ndarray::Array<T,2,2> fullWorkspaceT = ndarray::allocate(
        computeSize(maxConvolvedOrder), x.template getSize<0>()
    );
    ndarray::Array<T,2,-2> fullWorkspace = fullWorkspaceT.transpose();
    for (
        typename Impl::ItemVector::iterator itemIter = _impl->items.begin();
        itemIter != _impl->items.end();
        ++itemIter
    ) {
        itemIter->workspace = fullWorkspace[
            ndarray::view()(0, computeSize(itemIter->convolution->getRowOrder()))
        ];
    }
}

template class MultiShapeletMatrixBuilder<float>;
template class MultiShapeletMatrixBuilder<double>;

}} // namespace lsst::shapelet
