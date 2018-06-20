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
#include "pybind11/eigen.h"
#include "ndarray/pybind11.h"

#include "lsst/shapelet/GaussHermiteProjection.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

PYBIND11_PLUGIN(gaussHermiteProjection) {
    py::module::import("lsst.afw.geom");
    py::module mod("gaussHermiteProjection");

    py::class_<GaussHermiteProjection, std::shared_ptr<GaussHermiteProjection>> clsGaussHermiteProjection(
            mod, "GaussHermiteProjection");

    clsGaussHermiteProjection.def(py::init<int>(), "maxOrder"_a);

    clsGaussHermiteProjection.def("compute",
                                  (Eigen::MatrixXd (GaussHermiteProjection::*)(
                                          afw::geom::ellipses::Quadrupole const &, int inputOrder,
                                          afw::geom::ellipses::Quadrupole const &, int outputOrder) const) &
                                          GaussHermiteProjection::compute);
    clsGaussHermiteProjection.def("compute",
                                  (Eigen::MatrixXd (GaussHermiteProjection::*)(
                                          afw::geom::LinearTransform const &, int inputOrder,
                                          afw::geom::LinearTransform const &, int outputOrder) const) &
                                          GaussHermiteProjection::compute);
    clsGaussHermiteProjection.def(
            "compute", (Eigen::MatrixXd (GaussHermiteProjection::*)(Eigen::Matrix2d const &, int,
                                                                    Eigen::Matrix2d const &, int) const) &
                               GaussHermiteProjection::compute);

    clsGaussHermiteProjection.def("getMaxOrder", &GaussHermiteProjection::getMaxOrder);

    return mod.ptr();
}

}  // shapelet
}  // lsst
