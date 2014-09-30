// -*- LSST-C++ -*-
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

#include "lsst/shapelet/FunctorKeys.h"

namespace lsst { namespace shapelet {

ShapeletFunctionKey ShapeletFunctionKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc,
    std::string const & ellipseUnit,
    std::string const & coeffUnit,
    int order,
    BasisTypeEnum basisType
) {
    ShapeletFunctionKey result;
    result._order = order;
    result._basisType = basisType;
    result._ellipseKey = afw::table::EllipseKey::addFields(schema, name, doc, ellipseUnit);
    result._coefficientsKey =
        afw::table::ArrayKey<double>::addFields(schema, name, doc, coeffUnit, computeSize(order));
    return result;
}

ShapeletFunction ShapeletFunctionKey::get(afw::table::BaseRecord const & record) const {
    return ShapeletFunction(
        _order, _basisType,
        _ellipseKey.get(record),
        _coefficientsKey.get(record)
    );
}

void ShapeletFunctionKey::set(afw::table::BaseRecord & record, ShapeletFunction const & value) const {
    LSST_THROW_IF_NE(
        _order, value.getOrder(),
        pex::exceptions::InvalidParameterError,
        "FunctorKey order (%d) does not match ShapeletFunction order (%d)"
    );
    if (value.getBasisType() != _basisType) {
        ShapeletFunction v2(value);
        v2.changeBasisType(_basisType);
        _coefficientsKey.set(record, v2.getCoefficients());
    } else {
        _coefficientsKey.set(record, value.getCoefficients());
    }
    _ellipseKey.set(record, value.getEllipse());
}

}} // namespace lsst::shapelet
