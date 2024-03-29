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

#include "lsst/shapelet/GaussHermiteConvolution.h"
#include "lsst/shapelet/ShapeletFunction.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

void wrapGaussHermiteConvolution(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyGaussHermiteConvolution = py::class_<GaussHermiteConvolution, std::shared_ptr<GaussHermiteConvolution>>;
    wrappers.wrapType(PyGaussHermiteConvolution(wrappers.module, "GaussHermiteConvolution"), [](auto &mod, auto &cls) {
        cls.def(py::init<int, ShapeletFunction const &>(), "colOrder"_a, "psf"_a);

        cls.def("computeRowOrder", &GaussHermiteConvolution::computeRowOrder);
        cls.def("evaluate", &GaussHermiteConvolution::evaluate);
        cls.def("getColOrder", &GaussHermiteConvolution::getColOrder);
        cls.def("getRowOrder", &GaussHermiteConvolution::getRowOrder);
    });
}

}  // shapelet
}  // lsst
