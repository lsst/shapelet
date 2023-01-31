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

#include "lsst/shapelet/MultiShapeletBasis.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

void wrapMultiShapeletBasis(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyClass = py::class_<MultiShapeletBasisComponent, std::shared_ptr<MultiShapeletBasisComponent>>;

    static auto component =
            wrappers.wrapType(PyClass(wrappers.module, "MultiShapeletBasisComponent"), [](auto &mod, auto &cls) {
                cls.def(py::init<double, int, ndarray::Array<double const, 2, 2> const &>(),
                        "radius"_a, "order"_a, "matrix"_a);

                cls.def("getRadius", &MultiShapeletBasisComponent::getRadius);
                cls.def("getOrder", &MultiShapeletBasisComponent::getOrder);
                cls.def("getMatrix", &MultiShapeletBasisComponent::getMatrix);
            });

    using PyMultiShapeletBasis = py::class_<MultiShapeletBasis, std::shared_ptr<MultiShapeletBasis>>;

    wrappers.wrapType(PyMultiShapeletBasis(wrappers.module, "MultiShapeletBasis"), [](auto &mod, auto &cls) {
        cls.attr("Component") = component;

        cls.def(py::init<int>());
        cls.def(py::init<MultiShapeletBasis const &>());

        cls.def("getSize", &MultiShapeletBasis::getSize);
        cls.def("getComponentCount", &MultiShapeletBasis::getComponentCount);
        cls.def("addComponent", &MultiShapeletBasis::addComponent);
        cls.def("scale", &MultiShapeletBasis::scale);
        cls.def("normalize", &MultiShapeletBasis::normalize);
        cls.def("merge", &MultiShapeletBasis::merge);
        cls.def("makeFunction", &MultiShapeletBasis::makeFunction);
    });
}

}  // shapelet
}  // lsst
