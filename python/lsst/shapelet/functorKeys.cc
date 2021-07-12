/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/shapelet/FunctorKeys.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

PYBIND11_MODULE(functorKeys, mod) {
    py::module::import("lsst.afw.table");

    py::class_<ShapeletFunctionKey, std::shared_ptr<ShapeletFunctionKey>> clsShapeletFunctionKey(
        mod, "ShapeletFunctionKey");

    clsShapeletFunctionKey.def(py::init<>());
    clsShapeletFunctionKey.def(
            py::init<afw::table::EllipseKey const &, afw::table::ArrayKey<double> const &, BasisTypeEnum>(),
            "ellipse"_a, "coefficients"_a, "basisType"_a = HERMITE);
    clsShapeletFunctionKey.def(py::init<afw::table::SubSchema const &, BasisTypeEnum>(), "s"_a,
                               "basisType"_a = HERMITE);

    clsShapeletFunctionKey.def("__eq__", &ShapeletFunctionKey::operator==, py::is_operator());
    clsShapeletFunctionKey.def("__ne__", &ShapeletFunctionKey::operator!=, py::is_operator());

    clsShapeletFunctionKey.def_static("addFields", &ShapeletFunctionKey::addFields, "schema"_a, "name"_a,
                                      "doc"_a, "ellipseUnit"_a, "coeffUnit"_a, "order"_a,
                                      "basisType"_a = HERMITE);

    clsShapeletFunctionKey.def("get", &ShapeletFunctionKey::get);
    clsShapeletFunctionKey.def("set", &ShapeletFunctionKey::set);
    clsShapeletFunctionKey.def("isValid", &ShapeletFunctionKey::isValid);
    clsShapeletFunctionKey.def("getEllipse", &ShapeletFunctionKey::getEllipse);
    clsShapeletFunctionKey.def("getCoefficients", &ShapeletFunctionKey::getCoefficients);
    clsShapeletFunctionKey.def("getOrder", &ShapeletFunctionKey::getOrder);
    clsShapeletFunctionKey.def("getBasisType", &ShapeletFunctionKey::getBasisType);

    py::class_<MultiShapeletFunctionKey, std::shared_ptr<MultiShapeletFunctionKey>>
            clsMultiShapeletFunctionKey(mod, "MultiShapeletFunctionKey");

    clsMultiShapeletFunctionKey.def(py::init<>());
    clsMultiShapeletFunctionKey.def(py::init<afw::table::SubSchema const &, BasisTypeEnum>(), "s"_a,
                                    "basisType"_a = HERMITE);
    clsMultiShapeletFunctionKey.def(py::init<std::vector<std::shared_ptr<ShapeletFunctionKey>> const &>(),
                                    "components"_a);

    clsMultiShapeletFunctionKey.def_static("addFields", MultiShapeletFunctionKey::addFields, "schema"_a,
                                           "name"_a, "doc"_a, "ellipseUnit"_a, "coeffUnit"_a, "orders"_a,
                                           "basisType"_a = HERMITE);

    clsMultiShapeletFunctionKey.def("__eq__", &MultiShapeletFunctionKey::operator==, py::is_operator());
    clsMultiShapeletFunctionKey.def("__ne__", &MultiShapeletFunctionKey::operator!=, py::is_operator());
    clsMultiShapeletFunctionKey.def("__getitem__",
                                    [](MultiShapeletFunctionKey &self, int i) { return self[i]; });

    clsMultiShapeletFunctionKey.def("get", &MultiShapeletFunctionKey::get);
    clsMultiShapeletFunctionKey.def("set", &MultiShapeletFunctionKey::set);
    clsMultiShapeletFunctionKey.def("isValid", &MultiShapeletFunctionKey::isValid);
}

}  // shapelet
}  // lsst
