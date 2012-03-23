// -*- lsst-c++ -*-
%define meas_extensions_multiShapelet_DOCSTRING
"
Measurement algorithms using galaxy models built from multi-scale shapelets.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.daf.base", docstring=meas_extensions_multiShapelet_DOCSTRING) multiShapeletLib

%{
#include "lsst/pex/logging.h"
#include "lsst/meas/algorithms.h"
%}

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()
%import "lsst/meas/algorithms/algorithmsLib.i"
