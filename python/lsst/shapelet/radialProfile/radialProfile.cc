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

#include "lsst/shapelet/RadialProfile.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

void wrapRadialProfile(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyRadialProfil = py::class_<RadialProfile, std::shared_ptr<RadialProfile>> ;
    wrappers.wrapType(PyRadialProfil(wrappers.module, "RadialProfile"), [](auto &mod, auto &cls) {
        cls.def_static("get", &RadialProfile::get, py::return_value_policy::reference);
        cls.def("getName", &RadialProfile::getName);
        cls.def("_evaluate", (double (RadialProfile::*)(double) const) &RadialProfile::evaluate);
        cls.def("_evaluate", (ndarray::Array<double, 1, 1> (RadialProfile::*)(
                ndarray::Array<double const, 1, 1> const &) const) &
                RadialProfile::evaluate);
        cls.def("getMomentsRadiusFactor", &RadialProfile::getMomentsRadiusFactor);
        cls.def("getBasis", &RadialProfile::getBasis, "nComponents"_a, "maxRadius"_a = 0);
        cls.def("registerBasis", &RadialProfile::registerBasis);
    });
}

}  // shapelet
}  // lsst
