// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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
#ifndef MULTISHAPELET_HybridOptimizer_h_INCLUDED
#define MULTISHAPELET_HybridOptimizer_h_INCLUDED

#include "boost/noncopyable.hpp"

#include "ndarray/eigen.h"

#include "lsst/base.h"
#include "lsst/pex/config.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class Objective : private boost::noncopyable {
public:

    virtual void computeFunction(
        ndarray::Array<double const,1,1> const & parameters, 
        ndarray::Array<double,1,1> const & function
    ) = 0;

    virtual void computeDerivative(
        ndarray::Array<double const,1,1> const & parameters, 
        ndarray::Array<double const,1,1> const & function,
        ndarray::Array<double,2,-2> const & derivative
    ) = 0;

    virtual bool computePrior(
        ndarray::Array<double const,1,1> const & parameters, 
        double & value,
        ndarray::Array<double,1,1> const & gradient, 
        ndarray::Array<double,2,-2> const & hessian 
    ) { return false; }

    virtual ~Objective() {}

    int getFunctionSize() const { return _functionSize; }
    int getParameterSize() const { return _parameterSize; }

    Objective(int functionSize, int parameterSize) :
        _functionSize(functionSize), _parameterSize(parameterSize) {}

private:
    int const _functionSize;
    int const _parameterSize;
};

class HybridOptimizerControl {
public:
    LSST_CONTROL_FIELD(fTol, double, "stopping tolerance for function value");
    LSST_CONTROL_FIELD(gTol, double, "stopping tolerance for gradient value");
    LSST_CONTROL_FIELD(minStep, double, "minimum step size");
    LSST_CONTROL_FIELD(maxIter, int, "maximum number of iterations");
    LSST_CONTROL_FIELD(tau, double, "LM parameter (FIXME!)");
    LSST_CONTROL_FIELD(delta0, double, "BFGS parameter (FIXME!)");
    LSST_CONTROL_FIELD(useCholesky, bool, "whether to use Cholesky or Eigensystem factorization");

    HybridOptimizerControl() : 
        fTol(1E-8), gTol(1E-8), minStep(1E-8), maxIter(200), tau(1E-3), delta0(1.0), useCholesky(true) {}
};


/**
 *  @brief Hybrid Levenberg-Marquardt and BFGS Quasi-Newton optimizer.
 *
 *  This is a reimplementation and modification of the "Hybrid" option in Mike Jarvis'
 *  NLSolver class, which as of this writing can be found in the meas/algorithms package.
 *  That is in turn based on
 *  "Methods for Non-Linear Least Squares Problems", 2nd edition,
 *  April, 2004, K. Madsen, H.B. Nielsen, O. Tingleff,
 *  Informatics and Mathematical Modelling, Technical University of Denmark.
 *
 *  The main goal of the reimplementation is to expose the main loop to the user,
 *  adding them to inspect each step in detail.
 */

class HybridOptimizer {
public:

    typedef HybridOptimizerControl Control;

    enum MethodEnum { LM=0, BFGS=1 };

    enum StateFlags {
        WORKING          = 0x00, ///< No stopping condition met.
        SUCCESS_FTOL     = 0x01, ///< Function values are below tolerance (i.e. perfect fit).
        SUCCESS_GTOL     = 0x02, ///< Gradient values are below tolerance (at minimum).
        SUCCESS          = SUCCESS_FTOL | SUCCESS_GTOL, ///< Any success condition.
        FAILURE_MINSTEP  = 0x04, ///< Calculated step became too small.
        FAILURE_MINTRUST = 0x08, ///< Trust region became too small.
        FAILURE_MAXITER  = 0x10, ///< Too many iterations.
        FAILURE          = FAILURE_MINSTEP | FAILURE_MINTRUST | FAILURE_MAXITER ///< Any failure condition.
    };

    /**
     *  @brief Take a single step.
     *
     *  Return value is a bitwise OR of StateFlags.
     *
     *  This will never return FAILURE_MAXITER, as iterations are only counted within a call to run().
     *  It will attempt to take a step even if one of the SUCCESS conditions are met, but it will
     *  not do so if the FAILURE_MINSTEP or FAILURE_MINTRUST conditions are met, as this would
     *  usually just introduce round-off error to the problem.
     */
    int step();

    /**
     *  @brief Return the current state of the optimizer.
     *
     *  Return value is a bitwise OR of StateFlags.
     */
    int getState() const;    

    /**
     *  @brief Call step() in a loop until it succeeds or fails.
     *
     *  Return value is a bitwise OR of StateFlags.
     */
    int run();

    /// @brief Return which mode (Levenberg-Marquardt or Quasi-Newton BFGS) the optimizer is in.
    MethodEnum getMethod() const;

    /// @brief Return the squared norm of the function vector.
    double getChiSq() const;

    /**
     *  @brief Return the maximum absolute value (inf norm) of the function vector.
     *
     *  The SUCCESS_FTOL condition is met when this is less than the fTol control value.
     */
    double getFunctionInfNorm() const;

    /**
     *  @brief Return the maximum absolute value (inf norm) of the gradient vector.
     *
     *  The SUCCESS_GTOL condition is met when this is less than the gTol control value.
     */
    double getGradientInfNorm() const;

    /// @brief Return the LM 'mu' parameter that multiplies the diagonal term added to the Hessian.
    double getMu() const;

    /// @brief Return the BFGS trust region size.
    double getDelta() const;

    ndarray::Array<double const,1,1> getParameters() const;

    ndarray::Array<double const,1,1> getFunction() const;

    Control const & getControl() const;

    explicit HybridOptimizer(
        PTR(Objective) const & objective, 
        ndarray::Array<double const,1,1> const & parameters,
        Control const & ctrl = Control()
    );

    ~HybridOptimizer();

private:

    class Impl;

    PTR(Impl) _impl;
};

}}}} // namespace lsst::meas::extensions::multisShapelet

#endif // !MULTISHAPELET_HybridOptimizer_h_INCLUDED
