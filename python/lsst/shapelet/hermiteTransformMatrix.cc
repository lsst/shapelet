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

#include "lsst/shapelet/HermiteTransformMatrix.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

PYBIND11_PLUGIN(_hermiteTransformMatrix) {
    py::module mod("_hermiteTransformMatrix");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    py::class_<HermiteTransformMatrix, std::shared_ptr<HermiteTransformMatrix>> clsHermiteTransformMatrix(
            mod, "HermiteTransformMatrix");

    clsHermiteTransformMatrix.def(py::init<int>(), "order"_a);

    clsHermiteTransformMatrix.def(
            "compute", (Eigen::MatrixXd (HermiteTransformMatrix::*)(Eigen::Matrix2d const &) const) &
                               HermiteTransformMatrix::compute,
            "transform"_a);
    clsHermiteTransformMatrix.def(
            "compute",
            (Eigen::MatrixXd (HermiteTransformMatrix::*)(afw::geom::LinearTransform const &) const) &
                    HermiteTransformMatrix::compute,
            "transform"_a);
    clsHermiteTransformMatrix.def(
            "compute", (Eigen::MatrixXd (HermiteTransformMatrix::*)(Eigen::Matrix2d const &, int) const) &
                               HermiteTransformMatrix::compute,
            "transform"_a, "order"_a);
    clsHermiteTransformMatrix.def("compute", (Eigen::MatrixXd (HermiteTransformMatrix::*)(
                                                     afw::geom::LinearTransform const &, int) const) &
                                                     HermiteTransformMatrix::compute,
                                  "transform"_a, "order"_a);

    clsHermiteTransformMatrix.def("getCoefficientMatrix", &HermiteTransformMatrix::getCoefficientMatrix);
    clsHermiteTransformMatrix.def("getInverseCoefficientMatrix",
                                  &HermiteTransformMatrix::getInverseCoefficientMatrix);
    clsHermiteTransformMatrix.def("getOrder", &HermiteTransformMatrix::getOrder);

    return mod.ptr();
}
}
}  // lsst::shapelet
