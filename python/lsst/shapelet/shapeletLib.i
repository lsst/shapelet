// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
%define shapeletLib_DOCSTRING
"
Python interface to lsst::shapelet classes and functions
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.shapelet", docstring=shapeletLib_DOCSTRING) shapeletLib

%{
#   include "lsst/afw/geom.h"
#   include "lsst/afw/detection.h"
#   include "lsst/afw/image.h"
#   include "lsst/afw/cameraGeom.h"
#   include "lsst/pex/logging.h"
#   include "lsst/shapelet.h"
%}

%include "lsst/p_lsstSwig.i"
%include "std_vector.i"

%{
#include "lsst/afw/geom.h"
#define PY_ARRAY_UNIQUE_SYMBOL LSST_AFW_MATH_SHAPELETS_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
#include "lsst/afw/geom/ellipses/PyPixelRegion.h"
#include "lsst/afw/table.h"
%}

%init %{
    import_array();
%}
%pythoncode %{
    import numpy
%}

%include "ndarray.i"

%declareNumPyConverters(Eigen::VectorXd);
%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(Eigen::Matrix2d);
%declareNumPyConverters(ndarray::Array<double,1>);
%declareNumPyConverters(ndarray::Array<double const,1>);
%declareNumPyConverters(ndarray::Array<float,1,1>);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,1>);
%declareNumPyConverters(ndarray::Array<float,2,-1>);
%declareNumPyConverters(ndarray::Array<double,2,-1>);
%declareNumPyConverters(ndarray::Array<float,2,-2>);
%declareNumPyConverters(ndarray::Array<double,2,-2>);
%declareNumPyConverters(ndarray::Array<float const,1,1>);
%declareNumPyConverters(ndarray::Array<double const,1,1>);
%declareNumPyConverters(ndarray::Array<double const,2,-2>);
%declareNumPyConverters(ndarray::Array<double const,2,2>);
%declareNumPyConverters(ndarray::Array<double,3,-3>);

%feature(valuewrapper) lsst::shapelet::ShapeletFunction;
%feature(valuewrapper) lsst::shapelet::MultiShapeletFunction;

%template(MultiShapeletComponentList) std::vector<lsst::shapelet::ShapeletFunction>;

%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/table/tableLib.i"

%lsst_exceptions();

%copyctor lsst::shapelet::MultiShapeletBasis;
%ignore lsst::shapelet::ShapeletFunction::operator=;

%shared_ptr(lsst::shapelet::MultiShapeletBasis);

%include "lsst/shapelet/constants.h"
%include "lsst/shapelet/ConversionMatrix.h"
%include "lsst/shapelet/ShapeletFunction.h"
%include "lsst/shapelet/MultiShapeletFunction.h"
%include "lsst/shapelet/BasisEvaluator.h"
%include "lsst/shapelet/HermiteTransformMatrix.h"
%include "lsst/shapelet/GaussHermiteConvolution.h"
%include "lsst/shapelet/GaussHermiteProjection.h"
%include "lsst/shapelet/MultiShapeletBasis.h"

%include "lsst/shapelet/MatrixBuilder.h"
%define %instantiateMatrixBuilder(T, SUFFIX)
%template(MatrixBuilder##SUFFIX) lsst::shapelet::MatrixBuilder<T>;
%template(MatrixBuilderFactory##SUFFIX) lsst::shapelet::MatrixBuilderFactory<T>;
%template(MatrixBuilderWorkspace##SUFFIX) lsst::shapelet::MatrixBuilderWorkspace<T>;
%pythoncode %{
MatrixBuilder##SUFFIX.Factory = MatrixBuilderFactory##SUFFIX
MatrixBuilder##SUFFIX.Workspace = MatrixBuilderWorkspace##SUFFIX
MatrixBuilderFactory##SUFFIX.Workspace = MatrixBuilderWorkspace##SUFFIX
%}
%enddef

%instantiateMatrixBuilder(float, F)
%instantiateMatrixBuilder(double, D)

%extend lsst::shapelet::ShapeletFunction {

%pythoncode %{
def __reduce__(self):
    return (ShapeletFunction, (self.getOrder(), self.getBasisType(),
                               self.getEllipse(), self.getCoefficients()))
%}

}

%extend lsst::shapelet::MultiShapeletFunction {

%pythoncode %{
def __reduce__(self):
    return (MultiShapeletFunction, (list(self.getComponents()),))
%}

}

%extend lsst::shapelet::RadialProfile {
    %feature("shadow") evaluate %{
    def evaluate(self, r):
        if isinstance(r, numpy.ndarray):
            return $action(self, r.ravel()).reshape(r.shape)
        else:
            return $action(self, r)
    %}
}
%include "lsst/shapelet/RadialProfile.h"

%declareFunctorKey(ShapeletFunction, lsst::shapelet::ShapeletFunction)
%shared_ptr(lsst::shapelet::ShapeletFunctionKey)
%include "lsst/shapelet/FunctorKeys.h"
%useValueEquality(lsst::shapelet::ShapeletFunctionKey)
