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

#include <string>

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


MultiShapeletFunctionKey MultiShapeletFunctionKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc,
    std::string const & ellipseUnit,
    std::string const & coeffUnit,
    std::vector<int> const & orders,
    BasisTypeEnum basisType
) {
    MultiShapeletFunctionKey result;
    result._components.reserve(orders.size());
    for (std::size_t i = 0; i < orders.size(); ++i) {
        result._components.push_back(
            std::make_shared<ShapeletFunctionKey>(
                ShapeletFunctionKey::addFields(
                    schema,
                    schema[name][std::to_string(i)].getPrefix(),
                    doc,
                    ellipseUnit,
                    coeffUnit,
                    orders[i],
                    basisType
                )
            )
        );
    }
    return result;
}

MultiShapeletFunctionKey::MultiShapeletFunctionKey(
    afw::table::SubSchema const & s,
    BasisTypeEnum basisType
) {
    std::size_t i = 0;
    while (true) {
        try {
            std::shared_ptr<ShapeletFunctionKey> component = std::make_shared<ShapeletFunctionKey>(
                s[std::to_string(i)],
                basisType
            );
            _components.push_back(component);
            ++i;
        } catch (pex::exceptions::NotFoundError &) {
            break;
        }
    }
    if (i == 0) {
        throw LSST_EXCEPT(
            pex::exceptions::NotFoundError,
            (boost::format("No shapelet expansions found with prefix %s") % s.getPrefix()).str()
        );
    }
}

MultiShapeletFunction MultiShapeletFunctionKey::get(afw::table::BaseRecord const & record) const {
    MultiShapeletFunction result;
    result.getComponents().reserve(_components.size());
    for (std::size_t i = 0; i < _components.size(); ++i) {
        result.getComponents().push_back(_components[i]->get(record));
    }
    return result;
}

void MultiShapeletFunctionKey::set(
    afw::table::BaseRecord & record,
    MultiShapeletFunction const & value
) const {
    LSST_THROW_IF_NE(
        value.getComponents().size(), _components.size(),
        pex::exceptions::InvalidParameterError,
        "Number of components in value (%d) does not match number of components in FunctorKey (%d)"
    );
    for (std::size_t i = 0; i < _components.size(); ++i) {
        _components[i]->set(record, value.getComponents()[i]);
    }
}

bool MultiShapeletFunctionKey::operator==(MultiShapeletFunctionKey const & other) const {
    if (_components.size() != other._components.size()) {
        return false;
    }
    for (std::size_t i = 0; i < _components.size(); ++i) {
        if ((_components[i] != other._components[i]) && (*_components[i] != *other._components[i])) {
            return false;
        }
    }
    return true;
}

bool MultiShapeletFunctionKey::isValid() const {
    for (std::size_t i = 0; i < _components.size(); ++i) {
        if (!_components[i]->isValid()) return false;
    }
    return true;
}


}} // namespace lsst::shapelet
