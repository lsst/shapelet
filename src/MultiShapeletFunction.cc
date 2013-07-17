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

#include "lsst/shapelet/MultiShapeletFunction.h"
#include "lsst/shapelet/ConversionMatrix.h"
#include "lsst/pex/exceptions.h"
#include "ndarray/eigen.h"
#include <boost/format.hpp>

namespace lsst { namespace shapelet {

void MultiShapeletFunction::normalize(double value) {
    double factor = value / evaluate().integrate();
    for (ElementList::iterator i = _elements.begin(); i != _elements.end(); ++i) {
        i->getCoefficients().deep() *= factor;
    }
}

void MultiShapeletFunction::shiftInPlace(afw::geom::Extent2D const & offset) {
    for (ElementList::iterator i = _elements.begin(); i != _elements.end(); ++i) {
        i->shiftInPlace(offset);
    }    
}

void MultiShapeletFunction::transformInPlace(afw::geom::AffineTransform const & transform) {
    for (ElementList::iterator i = _elements.begin(); i != _elements.end(); ++i) {
        i->transformInPlace(transform);
    }    
}

MultiShapeletFunction MultiShapeletFunction::convolve(ShapeletFunction const & other) const {
    ElementList newElements;
    for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        newElements.push_back(i->convolve(other));
    }
    return MultiShapeletFunction(newElements);
}

MultiShapeletFunction MultiShapeletFunction::convolve(MultiShapeletFunction const & other) const {
    ElementList newElements;
    for (ElementList::const_iterator j = other.getElements().begin(); j != other.getElements().end(); ++j) {
        for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
            newElements.push_back(i->convolve(*j));
        }
    }
    return MultiShapeletFunction(newElements);
}

void MultiShapeletFunctionEvaluator::update(MultiShapeletFunction const & function) {
    _elements.clear();
    for (
        MultiShapeletFunction::ElementList::const_iterator i = function.getElements().begin(); 
        i != function.getElements().end();
        ++i
    ) {
        _elements.push_back(i->evaluate());
    }
}

double MultiShapeletFunctionEvaluator::operator()(afw::geom::Point2D const & point) const {
    double r = 0.0;
    for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        r += (*i)(point);
    }
    return r;
}

double MultiShapeletFunctionEvaluator::operator()(afw::geom::Extent2D const & point) const {
    double r = 0.0;
    for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        r += (*i)(point);
    }
    return r;
}

void MultiShapeletFunctionEvaluator::addToImage(
    ndarray::Array<double,2,1> const & array,
    afw::geom::Point2I const & xy0
) const {
    for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        i->addToImage(array, xy0);
    }
}

double MultiShapeletFunctionEvaluator::integrate() const {
    double r = 0.0;
    for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        r += i->integrate();
    }
    return r;
}

MultiShapeletFunctionEvaluator::MultiShapeletFunctionEvaluator(
    MultiShapeletFunction const & function
) {
    update(function);
}

afw::geom::ellipses::Ellipse MultiShapeletFunctionEvaluator::computeMoments() const {
    double q0 = 0.0;
    Eigen::Vector2d q1 = Eigen::Vector2d::Zero();
    Eigen::Matrix2d q2 = Eigen::Matrix2d::Zero();
    for (ElementList::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        i->_computeRawMoments(q0, q1, q2);
    }
    q1 /= q0;
    q2 /= q0;
    q2 -= q1 * q1.transpose();
    return afw::geom::ellipses::Ellipse(
        afw::geom::ellipses::Quadrupole(afw::geom::ellipses::Quadrupole::Matrix(q2), false),
        afw::geom::Point2D(q1)
    );
}

}} // namespace lsst::shapelet
