/*
 * This file is part of shapelet.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;
using lsst::cpputils::python::WrapperCollection;

namespace lsst {
namespace shapelet {

void wrapConstants(WrapperCollection &wrappers);
void wrapMatrixBuilder(WrapperCollection &wrappers);
void wrapGaussHermiteConvolution(WrapperCollection &wrappers);
void wrapFunctorKeys(WrapperCollection &wrappers);
void wrapShapeletFunction(WrapperCollection &wrappers);
void wrapGaussHermiteProjection(WrapperCollection &wrappers);
void wrapMultiShapeletFunction(WrapperCollection &wrappers);
void wrapHermiteTransformMatrix(WrapperCollection &wrappers);
void wrapMultiShapeletBasis(WrapperCollection &wrappers);
void wrapBasisEvaluator(WrapperCollection &wrappers);
void wrapRadialProfile(WrapperCollection &wrappers);

PYBIND11_MODULE(_shapeletLib, mod) {
    WrapperCollection wrappers(mod, "lsst.shapelet");

    wrappers.addInheritanceDependency("lsst.afw.table");

    wrappers.addSignatureDependency("lsst.geom");
    wrappers.addSignatureDependency("lsst.afw.geom");
    wrappers.addSignatureDependency("lsst.afw.image");

    wrapConstants(wrappers);
    wrapMatrixBuilder(wrappers);
    wrapGaussHermiteConvolution(wrappers);
    wrapFunctorKeys(wrappers);
    wrapShapeletFunction(wrappers);
    wrapGaussHermiteProjection(wrappers);
    wrapMultiShapeletFunction(wrappers);
    wrapHermiteTransformMatrix(wrappers);
    wrapMultiShapeletBasis(wrappers);
    wrapBasisEvaluator(wrappers);
    wrapRadialProfile(wrappers);
    wrappers.finish();
}

}  // namcespace shapelet
}  // namespace lsst