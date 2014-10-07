// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%module timeModels

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%{
#include "lsst/shapelet.h"
#include "lsst/afw/fits.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/table.h"
#include "lsst/pex/logging.h"
#define PY_ARRAY_UNIQUE_SYMBOL LSST_SHAPELET_TIMEMODELS_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%import "lsst/afw/table/tableLib.i"
%import "lsst/shapelet/shapeletLib.i"

%init %{
    import_array();
%}

%include "ndarray.i"

%declareNumPyConverters(ndarray::Array<double const,2>);

%{

#include "lsst/afw/geom/ellipses.h"

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
        ellipse.setParameterVector(i->asEigen());
        builder(matrix, ellipse);
        result += matrix.asEigen().norm();
    }
    return result;
}

%}

%inline %{

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

%}

%pythoncode %{

import numpy
import cPickle
import os
import sys
import argparse
import time
import resource
import lsst.shapelet.tractor

def main():
    parser = argparse.ArgumentParser(description="Benchmark shapelet-based galaxy model evaluation")
    parser.add_argument("-p", "--psf", help="Pickled PSF to use",
                        default="psf.p", type=str, metavar="FILE")
    parser.add_argument("-r", "--repeat", help="Number of times to repeat the test (loop in Python)",
                        default=100, type=int)
    parser.add_argument("-n", "--n-samples", help="Number of parameter points to evaluate (loop in C++)",
                        default=200, type=int)
    parser.add_argument("--n-disk-terms", help="Number of Gaussians in disk (exponential) component",
                        default=8, type=int)
    parser.add_argument("--n-bulge-terms", help="Number of Gaussians in bulge (de Vaucouleur) component",
                        default=8, type=int)
    parser.add_argument("--n-x-pixels", help="Number of pixels in x direction (footprint is a rectangle)",
                        default=20, type=int)
    parser.add_argument("--n-y-pixels", help="Number of pixels in y direction (footprint is a rectangle)",
                        default=20, type=int)
    parser.add_argument("--double-precision", help="Use double precision instead of single precision",
                        default=False, action="store_true")
    args = parser.parse_args()
    bulge = lsst.shapelet.tractor.loadBasis("luv", args.n_bulge_terms)
    disk = lsst.shapelet.tractor.loadBasis("lux", args.n_disk_terms)
    bulge.scale(0.6)
    disk.merge(bulge)
    with open(os.path.join(os.path.split(__file__)[0], args.psf), "r") as psfFile:
        psf = cPickle.load(psfFile)
    xMin = -args.n_x_pixels // 2
    yMin = -args.n_y_pixels // 2
    if args.double_precision:
        dtype = numpy.float64
        Builder = lsst.shapelet.MatrixBuilderD
        func = buildModelsD
    else:
        dtype = numpy.float32
        Builder = lsst.shapelet.MatrixBuilderF
        func = buildModelsF

    x, y = numpy.meshgrid(numpy.arange(xMin, xMin+args.n_x_pixels, dtype=dtype),
                          numpy.arange(yMin, yMin+args.n_y_pixels, dtype=dtype))
    builder = Builder(x.flatten(), y.flatten(), disk, psf)
    parameters = numpy.random.randn(args.n_samples, 5)
    cpuTime1 = time.clock()
    res1 = resource.getrusage(resource.RUSAGE_SELF)
    for n in range(args.repeat):
        func(x.size, disk.getSize(), builder, parameters, "SeparableConformalShearLogTraceRadius")
    cpuTime2 = time.clock()
    res2 = resource.getrusage(resource.RUSAGE_SELF)
    factor = args.n_samples * args.repeat
    print "Time per model evaluation (seconds): cpu=%g, user=%g, system=%g" % (
        (cpuTime2 - cpuTime1) / factor,
        (res2.ru_utime - res1.ru_utime) / factor,
        (res2.ru_stime - res1.ru_stime) / factor
        )

if __name__ == "__main__":
    main()
%}
