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

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/Cholesky"
#include "boost/make_shared.hpp"

#include "lsst/meas/extensions/multiShapelet/HybridOptimizer.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {
//
// Throughout this file, I have used the variable names from the original implementation
// and/or the formulae in that document, rather than names that adhere to the LSST naming
// conventions (which are follwed in the header file).  No amount of longer, more
// descriptive names would make this code understandable without a mathematical
// understanding of the algorithm, so I think it's more important to use variable names
// that can easily be mapped to variables in a more complete description.

class HybridOptimizer::Impl {
public:
    
    Impl(
        PTR(Objective) const & objective,
        ndarray::Array<double const,1,1> const & parameters,
        Control const & control
    );

    void step();

    void solve(Eigen::MatrixXd const & m);

    bool checkStep(double stepNorm, StateFlags bad) {
        if (!(stepNorm > ctrl.minStep * (x.norm() + ctrl.minStep))) {
            state |= bad;
            return false;
        }
        return true;
    }

    PTR(Objective) obj;
    HybridOptimizerControl ctrl;
    MethodEnum method;
    int state;
    int count;
    ndarray::EigenView<double,1,1> x;
    ndarray::EigenView<double,1,1> xNew;
    ndarray::EigenView<double,1,1> f;
    ndarray::EigenView<double,1,1> fNew;
    ndarray::EigenView<double,2,-2> J;
    ndarray::EigenView<double,2,-2> JNew;
    Eigen::VectorXd h;
    Eigen::VectorXd y;
    Eigen::VectorXd v;
    Eigen::MatrixXd A; // Hessian for LM method
    Eigen::MatrixXd B; // Hessian for BFGS method (Jarvis uses 'H')
    Eigen::VectorXd g;
    Eigen::VectorXd gNew;
    double normInfF;
    double normInfG;
    double Q;
    double QNew;
    double mu;
    double nu;
    double delta;
};

HybridOptimizer::Impl::Impl(
    PTR(Objective) const & objective,
    ndarray::Array<double const,1,1> const & parameters,
    Control const & control
) : obj(objective), ctrl(control), method(LM), state(WORKING), count(0),
    x(ndarray::copy(parameters)), xNew(ndarray::copy(parameters)),
    f(ndarray::allocate(objective->getFunctionSize())),
    fNew(ndarray::allocate(objective->getFunctionSize())),
    J(ndarray::allocate(objective->getFunctionSize(), objective->getParameterSize())),
    JNew(ndarray::allocate(objective->getFunctionSize(), objective->getParameterSize())),
    h(Eigen::VectorXd::Zero(objective->getParameterSize())),
    y(Eigen::VectorXd::Zero(objective->getParameterSize())),
    v(Eigen::VectorXd::Zero(objective->getParameterSize())),
    A(Eigen::VectorXd::Zero(objective->getParameterSize(), objective->getParameterSize())),
    B(Eigen::VectorXd::Identity(objective->getParameterSize(), objective->getParameterSize())),
    g(Eigen::VectorXd::Zero(objective->getFunctionSize())),
    gNew(Eigen::VectorXd::Zero(objective->getFunctionSize())),
    normInfF(0.0), normInfG(0.0), Q(0.0), QNew(0.0), mu(0.0), nu(2.0), delta(ctrl.delta0)
{
    fNew.setZero();
    obj->computeFunction(x.shallow(), f.shallow());
    normInfF = f.lpNorm<Eigen::Infinity>();
    Q = 0.5 * f.squaredNorm();
    JNew.setZero(); 
    obj->computeDerivative(x.shallow(), fNew.shallow(), JNew.shallow());
    f = fNew;
    J = JNew;
    A.selfadjointView<Eigen::Lower>().rankUpdate(J.adjoint());
    g = J.adjoint() * f;
    normInfG = g.lpNorm<Eigen::Infinity>();
    mu = ctrl.tau * A.diagonal().lpNorm<Eigen::Infinity>();
    A.diagonal().array() += mu;
}

void HybridOptimizer::Impl::step() {
    static double const sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());
    bool isBetter = false;
    bool shouldSwitchMethod = false;

    switch (method) {
    case LM:
        solve(A);
        break;
    case BFGS:
        solve(B);
        break;
    }
    
    double normH = h.norm();
    if (!checkStep(normH, FAILURE_MINSTEP)) return;
    if (method == BFGS && normH > delta) h *= delta / normH;
    
    xNew = x + h;
    fNew.setZero();
    obj->computeFunction(xNew.shallow(), fNew.shallow());
    QNew = 0.5 * fNew.squaredNorm();
    JNew.setZero();
    obj->computeDerivative(xNew.shallow(), fNew.shallow(), JNew.shallow());
    double normInfGNew = 0.0;
    if (method == BFGS || QNew < Q) {
        gNew = JNew.adjoint() * fNew;
        normInfGNew = gNew.lpNorm<Eigen::Infinity>();
    }

    if (method == BFGS) {
        isBetter = (QNew < Q) || (QNew <= (1.0 + sqrtEps) * Q && normInfGNew < normInfG);
        shouldSwitchMethod = (normInfGNew >= normInfG);
        if (QNew < Q) {
            double rho = (Q - QNew) / -(h.dot(g) - 0.5*(J*h).squaredNorm());
            if (rho > 0.75) {
                delta = std::max(delta, 3.0 * normH);
            } else if (rho < 0.25) {
                delta /= 2.0;
                if (!checkStep(delta, FAILURE_MINTRUST)) return;
            }
        } else {
            delta /= 2.0;
            if (!checkStep(delta, FAILURE_MINTRUST)) return;
        }
    } else { // method == LM
        if (QNew < Q) {
            isBetter = true;
            // replace g with g - mu*h, because we don't need it anymore
            g -= mu * h;
            double rho = (Q - QNew) / (-0.5 * h.dot(g));
            mu *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * rho - 1.0, 3));
            nu = 2.0;
            if (std::min(normInfGNew, Q - QNew) < 0.02 * QNew) {
                if (++count == 3) shouldSwitchMethod = true;
            } else {
                count = 0;
            }
            if (count != 3) {
                A.setZero();
                A.selfadjointView<Eigen::Lower>().rankUpdate(JNew.adjoint());
                A.diagonal().array() += mu;
            }
        } else {
            A.diagonal().array() += mu * (nu - 1.0);
            mu *= nu;
            nu *= 2.0;
            shouldSwitchMethod = (nu >= 32.0);
        }
    }

    y = JNew.adjoint() * (JNew * h) + (JNew - J).adjoint() * fNew;
    double hy = h.dot(y);
    if (hy > 0.0) {
        v = B.selfadjointView<Eigen::Lower>() * h;
        double hv = h.dot(v);
        B.selfadjointView<Eigen::Lower>().rankUpdate(v, -1.0 / hv);
        B.selfadjointView<Eigen::Lower>().rankUpdate(y, 1.0 / hy);
    }

    if (isBetter) {
        x = xNew;
        f = fNew;
        Q = QNew;
        J = JNew;
        g = gNew;
        normInfF = f.lpNorm<Eigen::Infinity>();
        normInfG = normInfGNew;
        if (!(normInfF > ctrl.fTol)) {
            state |= SUCCESS_FTOL;
        }
        if (!(normInfG > ctrl.gTol)) {
            state |= SUCCESS_GTOL;
        }
    }

    if (shouldSwitchMethod) {
        if (method == BFGS) { // switching from BFGS to LM
            A.setZero();
            A.selfadjointView<Eigen::Lower>().rankUpdate(J.adjoint());
            A.diagonal().array() += mu;
            method = LM;
        } else { // switching from LM to BFGS
            delta = std::max(1.5 * ctrl.minStep * (f.squaredNorm() + ctrl.minStep), 0.2 * normH);
            method = BFGS;
        }
    }
}


int HybridOptimizer::step() {
    _impl->step();
    return _impl->state;
}

int HybridOptimizer::getState() const {
    return _impl->state;
}

int HybridOptimizer::run() {
    for (int n = 0; n < _impl->ctrl.maxIter; ++n) {
        _impl->step();
        if (_impl->state) return _impl->state;
    }
    _impl->state |= FAILURE_MAXITER;
    return _impl->state;
}

HybridOptimizer::MethodEnum HybridOptimizer::getMethod() const { return _impl->method; }
double HybridOptimizer::getChiSq() const { return 2.0 * _impl->Q; }
double HybridOptimizer::getFunctionInfNorm() const { return _impl->normInfF; }
double HybridOptimizer::getGradientInfNorm() const { return _impl->normInfG; }
double HybridOptimizer::getMu() const { return _impl->mu; }
double HybridOptimizer::getDelta() const { return _impl->delta; }

ndarray::Array<double const,1,1> HybridOptimizer::getParameters() const { return _impl->x.shallow(); }
ndarray::Array<double const,1,1> HybridOptimizer::getFunction() const { return _impl->f.shallow(); }

HybridOptimizerControl const & HybridOptimizer::getControl() const { return _impl->ctrl; }

HybridOptimizer::HybridOptimizer(
    PTR(Objective) const & objective,
    ndarray::Array<double const,1,1> const & parameters,
    Control const & ctrl
) : _impl(boost::make_shared<Impl>(objective, parameters, ctrl))
{}

HybridOptimizer::~HybridOptimizer() {}

}}}} // namespace lsst::meas::extensions::multiShapelet
