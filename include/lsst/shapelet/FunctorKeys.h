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

#ifndef LSST_AFW_MATH_SHAPELETS_FunctorKeys_h_INCLUDED
#define LSST_AFW_MATH_SHAPELETS_FunctorKeys_h_INCLUDED

#include "lsst/afw/table/FunctorKey.h"
#include "lsst/afw/table/aggregates.h"
#include "lsst/afw/table/arrays.h"
#include "lsst/shapelet/ShapeletFunction.h"
#include "lsst/shapelet/MultiShapeletFunction.h"

namespace lsst { namespace shapelet {

/**
 *  Class that maps ShapeletFunction objects to fields in afw::table objects.
 *
 *  A ShapeletFunctionKey manages a set of fields with a common prefix and the following suffixes:
 *   - "x", "y", "xx", "yy", "xy": ellipse
 *   - "0", "1", "2", ...: coefficients.
 *
 *  As with all FunctorKeys, a ShapeletFunctorKey can be used to directly get or set objects
 *  on an afw::table::BaseRecord, just as with a true Key.
 */
class ShapeletFunctionKey : public afw::table::FunctorKey< ShapeletFunction > {
public:

    /**
     *  Add the necessary fields for saving a ShapeletFunction to a Schema.
     *
     *  @param[in,out] schema         Schema to add fields to.
     *  @param[in]     name           Name prefix for all fields.
     *  @param[in]     doc            String used as the documentation for the fields.
     *  @param[in]     ellipseUnit    String used as the unit for the ellipse ("<ellipseUnit>^2" will be
     *                                used for the Quadrupole moments).
     *  @param[in]     coeffUnit      String used as the unit for the coefficient vector
     *  @param[in]     order          Order of the ShapeletFunction to be saved.
     *  @param[in]     basisType      Type of shapelet basis (HERMITE or LAGUERRE) to be saved.
     *
     *  This method provides only basic exception safety - the schema may be (partially) modified even
     *  if an exception is thrown.
     */
    static ShapeletFunctionKey addFields(
        afw::table::Schema & schema,
        std::string const & name,
        std::string const & doc,
        std::string const & ellipseUnit,
        std::string const & coeffUnit,
        int order,
        BasisTypeEnum basisType=HERMITE
    );

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    ShapeletFunctionKey() : _ellipseKey(), _coefficientsKey(), _order(-1), _basisType(HERMITE) {}

    /// Construct from individual Keys/FunctorKeys
    ShapeletFunctionKey(
        afw::table::EllipseKey const & ellipse,
        afw::table::ArrayKey<double> const & coefficients,
        BasisTypeEnum basisType=HERMITE
    ) :
        _ellipseKey(ellipse),
        _coefficientsKey(coefficients),
        _order(computeOrder(coefficients.getSize())),
        _basisType(basisType)
    {}

    /**
     *  @brief Construct from a subschema, assuming the necesary subfields
     *
     *  If a schema has e.g. "a_xx", "a_0", etc. fields, this constructor allows you to construct
     *  a ShapeletFunctionKey via:
     *  @code
     *  ShapeletFunctionKey k(schema["a"]);
     *  @endcode
     */
    ShapeletFunctionKey(afw::table::SubSchema const & s, BasisTypeEnum basisType=HERMITE) :
        _ellipseKey(s), _coefficientsKey(s),
        _order(computeOrder(_coefficientsKey.getSize())),
        _basisType(basisType)
    {}

    /// Get a ShapeletFunction from the given record
    virtual ShapeletFunction get(afw::table::BaseRecord const & record) const;

    /// Set a ShapeletFunction in the given record
    virtual void set(afw::table::BaseRecord & record, ShapeletFunction const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Ixx, Iyy, Ixy Keys
    bool operator==(ShapeletFunctionKey const & other) const {
        return _order == other._order && _ellipseKey == other._ellipseKey
            && _coefficientsKey == other._coefficientsKey;
    }
    bool operator!=(ShapeletFunctionKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if all the constituent Keys are valid.
    bool isValid() const { return _ellipseKey.isValid() && _coefficientsKey.isValid(); }

    /// Return a FunctorKey that extracts just the Ellipse
    afw::table::EllipseKey const & getEllipse() const { return _ellipseKey; }

    /// Return a FunctorKey that extracts just the coefficients
    afw::table::ArrayKey<double> const & getCoefficients() const { return _coefficientsKey; }

    /// Return the shapelet order
    int getOrder() const { return _order; }

    /// Return the type of the shapelet basis
    BasisTypeEnum getBasisType() const { return _basisType; }

private:
    afw::table::EllipseKey _ellipseKey;
    afw::table::ArrayKey<double> _coefficientsKey;
    int _order;
    BasisTypeEnum _basisType;
};


/**
 *  Class that maps MultiShapeletFunction objects to fields in afw::table objects.
 *
 *  A MultiShapeletFunctionKey holds a sequnece of ShapeletFunctionKey, with an numerical prefix
 *  in front of each component.  A two-component MultiShapeletFunctionKey would thus be associated with
 *  the following keys:
 *   - "<prefix>_0_xx"
 *   - "<prefix>_0_yy"
 *   - "<prefix>_0_xy"
 *   - "<prefix>_0_x"
 *   - "<prefix>_0_y"
 *   - "<prefix>_1_xx"
 *   - "<prefix>_1_yy"
 *   - "<prefix>_1_xy"
 *   - "<prefix>_1_x"
 *   - "<prefix>_1_y"
 *
 *  As with all FunctorKeys, a MultiShapeletFunctorKey can be used to directly get or set objects
 *  on an afw::table::BaseRecord, just as with a true Key.
 */
class MultiShapeletFunctionKey : public afw::table::FunctorKey< MultiShapeletFunction > {
public:

    /**
     *  Add the necessary fields for saving a ShapeletFunction to a Schema.
     *
     *  @param[in,out] schema         Schema to add fields to.
     *  @param[in]     name           Name prefix for all fields.
     *  @param[in]     doc            String used as the documentation for the fields.
     *  @param[in]     ellipseUnit    String used as the unit for the ellipse ("<ellipseUnit>^2" will be
     *                                used for the Quadrupole moments).
     *  @param[in]     coeffUnit      String used as the unit for the coefficient vector
     *  @param[in]     orders         Vector of orders of the ShapeletFunctions to be saved.
     *  @param[in]     basisType      Type of shapelet basis (HERMITE or LAGUERRE) to be saved.
     *
     *  This method provides only basic exception safety - the schema may be (partially) modified even
     *  if an exception is thrown.
     */
    static MultiShapeletFunctionKey addFields(
        afw::table::Schema & schema,
        std::string const & name,
        std::string const & doc,
        std::string const & ellipseUnit,
        std::string const & coeffUnit,
        std::vector<int> const & orders,
        BasisTypeEnum basisType=HERMITE
    );

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    MultiShapeletFunctionKey() {}

    /// Construct from individual Keys/FunctorKeys
    explicit MultiShapeletFunctionKey(std::vector<std::shared_ptr<ShapeletFunctionKey>> const & components) :
        _components(components)
    {}

    /**
     *  @brief Construct from a subschema, assuming the necesary subfields
     *
     *  If a schema has e.g. "a_xx", "a_0", etc. fields, this constructor allows you to construct
     *  a ShapeletFunctionKey via:
     *  @code
     *  ShapeletFunctionKey k(schema["a"]);
     *  @endcode
     */
    MultiShapeletFunctionKey(afw::table::SubSchema const & s, BasisTypeEnum basisType=HERMITE);

    /// Get a MultiShapeletFunction from the given record
    virtual MultiShapeletFunction get(afw::table::BaseRecord const & record) const;

    /// Set a MultiShapeletFunction in the given record
    virtual void set(afw::table::BaseRecord & record, MultiShapeletFunction const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Ixx, Iyy, Ixy Keys
    bool operator==(MultiShapeletFunctionKey const & other) const;
    bool operator!=(MultiShapeletFunctionKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if all the constituent Keys are valid.
    bool isValid() const;

    /// Return a FunctorKey to the nth component.
    std::shared_ptr<ShapeletFunctionKey> operator[](int n) { return _components[n]; }

    /// Return a FunctorKey to the nth component.
    std::shared_ptr<ShapeletFunctionKey const> operator[](int n) const { return _components[n]; }

private:
    std::vector<std::shared_ptr<ShapeletFunctionKey>> _components;
};

}} // namespace lsst::shapelet

#endif // !defined(LSST_AFW_MATH_SHAPELETS_FunctorKeys_h_INCLUDED)
