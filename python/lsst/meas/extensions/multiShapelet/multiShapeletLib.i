// -*- lsst-c++ -*-
%define meas_extensions_multiShapelet_DOCSTRING
"
Measurement algorithms using galaxy models built from multi-scale shapelets.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.extensions.multiShapelet", 
        docstring=meas_extensions_multiShapelet_DOCSTRING) multiShapeletLib

%{
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection.h"
#include "lsst/pex/logging.h"
#include "lsst/meas/algorithms.h"
#include "lsst/meas/extensions/multiShapelet.h"
%}

%include "lsst/p_lsstSwig.i"

%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_AFW_MATH_SHAPELETS_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "ndarray.i"

%lsst_exceptions()
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/meas/algorithms/algorithmsLib.i"

%declareNumPyConverters(ndarray::Array<double const,1,1>);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,-1>);
%declareNumPyConverters(Eigen::Matrix<double,5,Eigen::Dynamic>);

%include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

%shared_ptr(lsst::meas::extensions::multiShapelet::Objective);
%include "lsst/meas/extensions/multiShapelet/HybridOptimizer.h"
