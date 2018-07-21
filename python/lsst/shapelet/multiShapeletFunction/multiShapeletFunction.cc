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

#include "ndarray/pybind11.h"

#include "lsst/shapelet/MultiShapeletFunction.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

namespace {

template <typename PyClass>
void declareMultiShapeletFunctionMembers(PyClass &cls) {
    using Class = MultiShapeletFunction;

    cls.def(py::init<>());
    cls.def(py::init<Class const &>());
    cls.def(py::init<typename Class::ComponentList const &>());
    cls.def(py::init<ShapeletFunction const &>());

    cls.def("getComponents", [](Class const &self) {
        auto &components = self.getComponents();

        py::tuple t(components.size());
        for (size_t i = 0; i < components.size(); ++i) {
            t[i] = py::cast(components[i]);
        }

        return t;
    });
    cls.def("addComponent",
            [](Class &self, typename Class::Component const &c) { self.getComponents().push_back(c); });
    cls.def("normalize", &Class::normalize, "value"_a = 1.0);
    cls.def("shiftInPlace", &Class::shiftInPlace);
    cls.def("transformInPlace", &Class::transformInPlace);
    cls.def("convolve", (Class (Class::*)(ShapeletFunction const &) const) & Class::convolve);
    cls.def("convolve", (Class (Class::*)(Class const &) const) & Class::convolve);
    cls.def("evaluate", &Class::evaluate);
}

template <typename PyClass>
void declareMultiShapeletFunctionEvaluatorMembers(PyClass &cls) {
    using Class = MultiShapeletFunctionEvaluator;

    cls.def(py::init<MultiShapeletFunction const &>());

    cls.def("__call__", (double (Class::*)(double, double) const) & Class::operator());
    cls.def("__call__", (double (Class::*)(afw::geom::Point2D const &) const) & Class::operator());
    cls.def("__call__", (double (Class::*)(afw::geom::Extent2D const &) const) & Class::operator());
    cls.def("__call__",
            (ndarray::Array<double, 1, 1> (Class::*)(ndarray::Array<double const, 1> const &,
                                                     ndarray::Array<double const, 1> const &) const) &
                    Class::operator());

    cls.def("addToImage",
            (void (Class::*)(ndarray::Array<double, 2, 1> const &, afw::geom::Point2I const &) const) &
                    Class::addToImage,
            "array"_a, "xy0"_a = afw::geom::Point2I());
    cls.def("addToImage", (void (Class::*)(afw::image::Image<double> &) const) & Class::addToImage,
            "image"_a);

    cls.def("integrate", &Class::integrate);
    cls.def("computeMoments", &Class::computeMoments);
    cls.def("update", &Class::update);
}

}  // <anonymous>

PYBIND11_MODULE(multiShapeletFunction, mod) {
    py::module::import("lsst.afw.geom");
    py::module::import("lsst.afw.image");

        py::class_<MultiShapeletFunction, std::shared_ptr<MultiShapeletFunction>> clsMultiShapeletFunction(
            mod, "MultiShapeletFunction");
    py::class_<MultiShapeletFunctionEvaluator, std::shared_ptr<MultiShapeletFunctionEvaluator>>
            clsMultiShapeletFunctionEvaluator(mod, "MultiShapeletFunctionEvaluator");

    declareMultiShapeletFunctionMembers(clsMultiShapeletFunction);
    declareMultiShapeletFunctionEvaluatorMembers(clsMultiShapeletFunctionEvaluator);
}

}  // shapelet
}  // lsst
