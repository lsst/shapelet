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

#include "ndarray/pybind11.h"

#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/shapelet/MultiShapeletBasis.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

namespace {

template <typename T>
py::class_<MatrixBuilder<T>, std::shared_ptr<MatrixBuilder<T>>> declareMatrixBuilder(
        py::module &mod, std::string const &suffix) {
    using Class = MatrixBuilder<T>;

    py::class_<Class, std::shared_ptr<Class>> cls(mod, ("MatrixBuilder" + suffix).c_str());

    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &, int>(),
            "x"_a, "y"_a, "order"_a);
    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &, int,
                     ShapeletFunction const &>(),
            "x"_a, "y"_a, "order"_a, "psf"_a);
    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &,
                     MultiShapeletBasis const &>(),
            "x"_a, "y"_a, "basis"_a);
    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &,
                     MultiShapeletBasis const &, MultiShapeletFunction const &>(),
            "x"_a, "y"_a, "basis"_a, "psf"_a);

    cls.def("getDataSize", &Class::getDataSize);
    cls.def("getBasisSize", &Class::getBasisSize);
    cls.def("allocateOutput", &Class::allocateOutput);

    cls.def("__call__",
            (void (Class::*)(ndarray::Array<T, 2, -1> const &, afw::geom::ellipses::Ellipse const &) const) &
                    Class::operator());
    cls.def("__call__", (ndarray::Array<T, 2, -2> (Class::*)(afw::geom::ellipses::Ellipse const &) const) &
                                Class::operator());

    return cls;
}

template <typename T>
py::class_<MatrixBuilderWorkspace<T>, std::shared_ptr<MatrixBuilderWorkspace<T>>>
declareMatrixBuilderWorkspace(py::module &mod, std::string const &suffix) {
    using Class = MatrixBuilderWorkspace<T>;

    py::class_<Class, std::shared_ptr<Class>> cls(mod, ("MatrixBuilderWorkspace" + suffix).c_str());

    cls.def(py::init<int>(), "size"_a);
    cls.def(py::init<Class const &>(), "other"_a);

    cls.def("getRemaining", &Class::getRemaining);

    return cls;
}

template <typename T>
py::class_<MatrixBuilderFactory<T>, std::shared_ptr<MatrixBuilderFactory<T>>> declareMatrixBuilderFactory(
        py::module &mod, std::string const &suffix) {
    using Class = MatrixBuilderFactory<T>;

    py::class_<Class, std::shared_ptr<Class>> cls(mod, ("MatrixBuilderFactory" + suffix).c_str());

    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &, int>(),
            "x"_a, "y"_a, "order"_a);
    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &, int,
                     ShapeletFunction const &>(),
            "x"_a, "y"_a, "order"_a, "psf"_a);
    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &,
                     MultiShapeletBasis const &>(),
            "x"_a, "y"_a, "basis"_a);
    cls.def(py::init<ndarray::Array<T const, 1, 1> const &, ndarray::Array<T const, 1, 1> const &,
                     MultiShapeletBasis const &, MultiShapeletFunction const &>(),
            "x"_a, "y"_a, "basis"_a, "psf"_a);

    cls.def("__call__", (MatrixBuilder<T> (Class::*)() const) & Class::operator());
    cls.def("__call__", (MatrixBuilder<T> (Class::*)(typename Class::Workspace &) const) & Class::operator(),
            "workspace"_a);

    cls.def("getDataSize", &Class::getDataSize);
    cls.def("getBasisSize", &Class::getBasisSize);
    cls.def("computeWorkspace", &Class::computeWorkspace);

    return cls;
}

template <typename T>
void declareMatrixBuilderTemplates(py::module &mod, std::string const &suffix) {
    auto clsMatrixBuilder = declareMatrixBuilder<T>(mod, suffix);
    auto clsMatrixBuilderWorkspace = declareMatrixBuilderWorkspace<T>(mod, suffix);
    auto clsMatrixBuilderFactory = declareMatrixBuilderFactory<T>(mod, suffix);

    clsMatrixBuilder.attr("Workspace") = clsMatrixBuilderWorkspace;
    clsMatrixBuilder.attr("Factory") = clsMatrixBuilderFactory;

    clsMatrixBuilderFactory.attr("Workspace") = clsMatrixBuilderWorkspace;
    clsMatrixBuilderFactory.attr("Builder") = clsMatrixBuilder;
}

}  // <anonymous>

PYBIND11_MODULE(matrixBuilder, mod) {
    py::module::import("lsst.afw.geom");
        declareMatrixBuilderTemplates<float>(mod, "F");
    declareMatrixBuilderTemplates<double>(mod, "D");
}

}  // shapelet
}  // lsst
