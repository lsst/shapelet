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
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include "ndarray/pybind11.h"

#include "lsst/shapelet.h"
#include "lsst/afw/fits.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/table.h"
#include "lsst/afw/geom/ellipses.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace shapelet {

namespace {
template <typename T>
double buildModelsImpl(
    int nPixels, int basisSize,
    lsst::shapelet::MatrixBuilder<T> const & builder,
    ndarray::Array<double const,2> const & parameters,
    std::string const & ellipseType
) {
    lsst::afw::geom::ellipses::Ellipse ellipse(
        lsst::afw::geom::ellipses::BaseCore::make(ellipseType),
        lsst::afw::geom::Point2D()
    );
    ndarray::Array<T,2,2> matrixT = ndarray::allocate(basisSize, nPixels);
    ndarray::Array<T,2,-2> matrix = matrixT.transpose();
    double result = 0.0;
    for (ndarray::Array<double const,2>::Iterator i = parameters.begin(); i != parameters.end(); ++i) {
        ellipse.setParameterVector(ndarray::asEigenMatrix(*i));
        builder(matrix, ellipse);
        result += ndarray::asEigenMatrix(matrix).norm();
    }
    return result;
}
    
void buildModelsF(
    int nPixels, int basisSize,
    lsst::shapelet::MatrixBuilder<float> const & builder,
    ndarray::Array<double const,2> const & parameters,
    std::string const & ellipseType
) {
    buildModelsImpl(nPixels, basisSize, builder, parameters, ellipseType);
}

void buildModelsD(
    int nPixels, int basisSize,
    lsst::shapelet::MatrixBuilder<double> const & builder,
    ndarray::Array<double const,2> const & parameters,
    std::string const & ellipseType
) {
    buildModelsImpl(nPixels, basisSize, builder, parameters, ellipseType);
}
}

PYBIND11_MODULE(_timeModels, mod) {
    mod.def("buildModelsF", buildModelsF);
    mod.def("buildModelsD", buildModelsD);
}

}} // lsst::shapelet
