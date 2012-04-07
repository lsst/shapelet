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
 
%define testLib_DOCSTRING
"
Various swigged-up C++ classes for testing
"
%enddef

%feature("autodoc", "1");
%module(package="testLib", docstring=testLib_DOCSTRING) testLib

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
#define PY_ARRAY_UNIQUE_SYMBOL LSST_TESTLIB_NUMPY_ARRAY_API
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

%import "lsst/meas/extensions/multiShapelet/multiShapeletLib.i"

%declareNumPyConverters(ndarray::Array<double const,1,1>);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,-2>);

%shared_ptr(RosenbrockObjective);

%inline %{

    class RosenbrockObjective : public lsst::meas::extensions::multiShapelet::Objective {
    public:

        explicit RosenbrockObjective(double lambda_) : 
            lsst::meas::extensions::multiShapelet::Objective(3, 2),
            _lambda(lambda_)
        {}

        virtual void computeFunction(
            ndarray::Array<double const,1,1> const & parameters, 
            ndarray::Array<double,1,1> const & function
        ) {
            function[0] = 10.0 * (parameters[1] - parameters[0] * parameters[0]);
            function[1] = 1.0 - parameters[0];
            function[2] = _lambda;
        }

        virtual void computeDerivative(
            ndarray::Array<double const,1,1> const & parameters, 
            ndarray::Array<double const,1,1> const & function,
            ndarray::Array<double,2,-2> const & derivative
        ) {
            derivative.deep() = 0.0;
            derivative[0][0] = -20.0 * parameters[0];
            derivative[0][1] = 10.0;
            derivative[1][0] = -1.0;
        }

    private:
        double _lambda;
    };

%}
