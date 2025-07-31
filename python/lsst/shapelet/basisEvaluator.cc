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

#include "lsst/shapelet/BasisEvaluator.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

void wrapBasisEvaluator(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyBasisEvaluator = py::classh<BasisEvaluator>;

    wrappers.wrapType(PyBasisEvaluator(wrappers.module, "BasisEvaluator"), [](auto &mod, auto &cls) {
        /* Constructors */
        cls.def(py::init<int, BasisTypeEnum>(), "order"_a, "basisType"_a);

        /* Members */
        cls.def("getOrder", &BasisEvaluator::getOrder);
        cls.def("getBasisType", &BasisEvaluator::getBasisType);

        /* fillEvaluation has default constructed Array1d objects as keyword
         * arguments.
         * Unfortunately, for some reason array.isEmpty() == false even with 0
         * elements,
         * which leads to a segfault.
         * Thus we must delegate through lambdas instead. */
        cls.def("fillEvaluation", [](BasisEvaluator &self, Array1d const &array, double x,
                                                   double y) { return self.fillEvaluation(array, x, y); },
                              "array"_a, "x"_a, "y"_a);
        cls.def(
                "fillEvaluation",
                [](BasisEvaluator &self, Array1d const &array, double x, double y, Array1d const &dx,
                   Array1d const &dy) { return self.fillEvaluation(array, x, y, dx, dy); },
                "array"_a, "x"_a, "y"_a, "dx"_a, "dy"_a);
        cls.def("fillEvaluation",
                              [](BasisEvaluator &self, Array1d const &array, geom::Point2D const &point) {
                                  return self.fillEvaluation(array, point);
                              },
                              "array"_a, "point"_a);
        cls.def(
                "fillEvaluation",
                [](BasisEvaluator &self, Array1d const &array, geom::Point2D const &point, Array1d const &dx,
                   Array1d const &dy) { return self.fillEvaluation(array, point, dx, dy); },
                "array"_a, "point"_a, "dx"_a, "dy"_a);
        cls.def("fillEvaluation",
                              [](BasisEvaluator &self, Array1d const &array, geom::Extent2D const &extent) {
                                  return self.fillEvaluation(array, extent);
                              },
                              "array"_a, "extent"_a);
        cls.def(
                "fillEvaluation",
                [](BasisEvaluator &self, Array1d const &array, geom::Extent2D const &extent,
                   Array1d const &dx, Array1d const &dy) { return self.fillEvaluation(array, extent, dx, dy); },
                "array"_a, "extent"_a, "dx"_a, "dy"_a);
        cls.def("fillIntegration", &BasisEvaluator::fillIntegration, "array"_a, "xMoment"_a = 0,
                              "yMoment"_a = 0);
    });
}

}  // shapelet
}  // lsst
