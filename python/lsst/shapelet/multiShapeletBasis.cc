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

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/shapelet/MultiShapeletBasis.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

PYBIND11_PLUGIN(multiShapeletBasis) {
    py::module mod("multiShapeletBasis");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    py::class_<MultiShapeletBasisComponent, std::shared_ptr<MultiShapeletBasisComponent>>
            clsMultiShapeletBasisComponent(mod, "MultiShapeletBasisComponent");

    clsMultiShapeletBasisComponent.def(py::init<double, int, ndarray::Array<double const, 2, 2> const &>(),
                                       "radius"_a, "order"_a, "matrix"_a);

    clsMultiShapeletBasisComponent.def("getRadius", &MultiShapeletBasisComponent::getRadius);
    clsMultiShapeletBasisComponent.def("getOrder", &MultiShapeletBasisComponent::getOrder);
    clsMultiShapeletBasisComponent.def("getMatrix", &MultiShapeletBasisComponent::getMatrix);

    py::class_<MultiShapeletBasis, std::shared_ptr<MultiShapeletBasis>> clsMultiShapeletBasis(
            mod, "MultiShapeletBasis");

    clsMultiShapeletBasis.attr("Component") = clsMultiShapeletBasisComponent;

    clsMultiShapeletBasis.def(py::init<int>());
    clsMultiShapeletBasis.def(py::init<MultiShapeletBasis const &>());

    clsMultiShapeletBasis.def("getSize", &MultiShapeletBasis::getSize);
    clsMultiShapeletBasis.def("getComponentCount", &MultiShapeletBasis::getComponentCount);
    clsMultiShapeletBasis.def("addComponent", &MultiShapeletBasis::addComponent);
    clsMultiShapeletBasis.def("scale", &MultiShapeletBasis::scale);
    clsMultiShapeletBasis.def("normalize", &MultiShapeletBasis::normalize);
    clsMultiShapeletBasis.def("merge", &MultiShapeletBasis::merge);
    clsMultiShapeletBasis.def("makeFunction", &MultiShapeletBasis::makeFunction);

    return mod.ptr();
}

}  // shapelet
}  // lsst
