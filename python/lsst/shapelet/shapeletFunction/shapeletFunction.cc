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
#include "lsst/cpputils/python.h"

#include "ndarray/pybind11.h"

#include "lsst/shapelet/ShapeletFunction.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

void wrapShapeletFunction(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyShapeletFunction = py::class_<ShapeletFunction>;

    wrappers.wrapType(PyShapeletFunction(wrappers.module, "ShapeletFunction"), [](auto &mod, auto &cls) {
        cls.def_readonly_static("FLUX_FACTOR", &ShapeletFunction::FLUX_FACTOR);

        cls.def(py::init<int, BasisTypeEnum>(), "order"_a, "basisType"_a);
        cls.def(py::init<int, BasisTypeEnum, ndarray::Array<double, 1, 1> const &>(), "order"_a,
                                "basisType"_a, "coefficients"_a);
        cls.def(py::init<int, BasisTypeEnum, double, geom::Point2D const &>(), "order"_a,
                                "basisType"_a, "radius"_a, "center"_a = geom::Point2D());
        cls.def(py::init<int, BasisTypeEnum, double, geom::Point2D const &,
                                        ndarray::Array<double, 1, 1> const &>(),
                                "order"_a, "basisType"_a, "radius"_a, "center"_a, "coefficients"_a);
        cls.def(py::init<int, BasisTypeEnum, afw::geom::ellipses::Ellipse const &>(), "order"_a,
                                "basisType"_a, "ellipse"_a);
        cls.def(py::init<int, BasisTypeEnum, afw::geom::ellipses::Ellipse const &,
                                        ndarray::Array<double const, 1, 1> const &>(),
                                "order"_a, "basisType"_a, "ellipse"_a, "coefficients"_a);
        cls.def(py::init<ShapeletFunction>());

        cls.def("getOrder", &ShapeletFunction::getOrder);
        cls.def("getEllipse", (afw::geom::ellipses::Ellipse &(ShapeletFunction::*)()) &
                                        ShapeletFunction::getEllipse,
                                py::return_value_policy::reference_internal);
        cls.def("setEllipse", &ShapeletFunction::setEllipse);
        cls.def("getBasisType", &ShapeletFunction::getBasisType);
        cls.def("changeBasisType", &ShapeletFunction::changeBasisType);
        cls.def("normalize", &ShapeletFunction::normalize, "value"_a = 1.0);
        cls.def("getCoefficients", (ndarray::Array<double, 1, 1> const (ShapeletFunction::*)()) &
                ShapeletFunction::getCoefficients);
        cls.def("convolve", &ShapeletFunction::convolve);
        cls.def("evaluate", &ShapeletFunction::evaluate);
        cls.def("shiftInPlace", &ShapeletFunction::shiftInPlace);
        cls.def("transformInPlace", &ShapeletFunction::transformInPlace);
    });

    using PyShapeletFunctionEvaluator = py::class_<ShapeletFunctionEvaluator>;

    wrappers.wrapType(
            PyShapeletFunctionEvaluator(wrappers.module, "ShapeletFunctionEvaluator"), [](auto &mod, auto &cls) {
                cls.def(py::init<ShapeletFunction const &>(), "function"_a);

                cls.def("__call__",
                        (double (ShapeletFunctionEvaluator::*)(double, double) const) &
                                ShapeletFunctionEvaluator::operator());
                cls.def(
                        "__call__", (double (ShapeletFunctionEvaluator::*)(geom::Point2D const &) const) &
                                ShapeletFunctionEvaluator::operator());
                cls.def(
                        "__call__", (double (ShapeletFunctionEvaluator::*)(geom::Extent2D const &) const) &
                                ShapeletFunctionEvaluator::operator());
                cls.def("__call__", (ndarray::Array<double, 1, 1> (ShapeletFunctionEvaluator::*)(
                        ndarray::Array<double const, 1> const &,
                        ndarray::Array<double const, 1> const &) const) &
                        ShapeletFunctionEvaluator::operator());

                cls.def(
                        "addToImage", (void (ShapeletFunctionEvaluator::*)(ndarray::Array<double, 2, 1> const &,
                                                                           geom::Point2I const &) const) &
                                ShapeletFunctionEvaluator::addToImage,
                        "array"_a, "xy0"_a = geom::Point2D());
                cls.def(
                        "addToImage", (void (ShapeletFunctionEvaluator::*)(afw::image::Image<double> &) const) &
                                ShapeletFunctionEvaluator::addToImage,
                        "image"_a);
                cls.def("integrate", &ShapeletFunctionEvaluator::integrate);
                cls.def("computeMoments", &ShapeletFunctionEvaluator::computeMoments);
                cls.def("update", &ShapeletFunctionEvaluator::update);
            });
}

}  // shapelet
}  // lsst
