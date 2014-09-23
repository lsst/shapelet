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

#include "ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

#include "lsst/shapelet/MultiShapeletBasis.h"

namespace lsst { namespace shapelet {

MultiShapeletBasisComponent::MultiShapeletBasisComponent(
    double radius,
    int order,
    ndarray::Array<double const,2,2> const & matrix
) : _radius(radius), _order(order), _matrix(ndarray::copy(matrix)) {
    if (_matrix.getSize<0>() != computeSize(order)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
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
            pex::exceptions::LengthError,
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
        result.getComponents().push_back(ShapeletFunction(i->getOrder(), HERMITE, ellipse));
        result.getComponents().back().getEllipse().getCore().scale(i->getRadius());
        result.getComponents().back().getCoefficients().asEigen()
            = i->getMatrix().asEigen() * coefficients.asEigen();
    }
    return result;
}

}} // namespace lsst::shapelet
