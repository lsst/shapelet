from __future__ import absolute_import, division, print_function

import numpy
import pickle
import os
import argparse
import time
import resource
import lsst.shapelet.tractor

from _timeModels import buildModelsD, buildModelsF


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
        psf = pickle.load(psfFile)
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
    print("Time per model evaluation (seconds): cpu=%g, user=%g, system=%g" % (
        (cpuTime2 - cpuTime1) / factor,
        (res2.ru_utime - res1.ru_utime) / factor,
        (res2.ru_stime - res1.ru_stime) / factor
    ))


if __name__ == "__main__":
    main()
