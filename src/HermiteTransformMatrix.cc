// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#include "lsst/shapelet/HermiteTransformMatrix.h"

namespace lsst { namespace shapelet {

namespace {

class Binomial {
public:

    Binomial(int const nMax);

    void reset(int const n, double a, double b);

    double const operator[](int const k) const { return _workspace[k] * _coefficients[k]; }

    int const getOrder() const { return _coefficients.size()-1; }

private:
    Eigen::VectorXd _coefficients;
    Eigen::VectorXd _workspace;
};

Binomial::Binomial(int const nMax) : _coefficients(nMax+1), _workspace(nMax+1) {}

void Binomial::reset(int const n, double a, double b) {
    assert(n < _coefficients.size());
    _coefficients[0] = _coefficients[n] = 1.0;
    int const mid = n/2;
    for (int k = 1; k <= mid; ++k) {
        _coefficients[k] = _coefficients[k-1] * (n - k + 1.0) / k;
    }
    for (int k = mid+1; k < n; ++k) {
        _coefficients[k] = _coefficients[n-k];
    }
    double v = 1;
    for (int k = 0; k <= n; ++k) {
        _workspace[k] = v;
        v *= b;
    }
    v = 1;
    for (int nk = n; nk >= 0; --nk) {
        _workspace[nk] *= v;
        v *= a;
    }
}

} // anonymous

HermiteTransformMatrix::HermiteTransformMatrix(int order) :
    _order(order),
    _coeffFwd(Eigen::MatrixXd::Zero(order+1, order+1)),
    _coeffInv(Eigen::MatrixXd::Identity(order+1, order+1))
{
    _coeffFwd(0, 0) = BASIS_NORMALIZATION;
    if (_order >= 1) {
        _coeffFwd(1, 1) = _coeffFwd(0, 0) * M_SQRT2;
    }
    for (int n = 2; n <= _order; ++n) {
        _coeffFwd(n, 0) = -_coeffFwd(n-2, 0) * rationalSqrt(n - 1, n);
        for (int m = (n % 2) ? 1:2; m <= n; m += 2) {
            _coeffFwd(n, m)
                = _coeffFwd(n-1, m-1) * rationalSqrt(2, n)
                - _coeffFwd(n-2, m) * rationalSqrt(n - 1, n);
        }
    }
    _coeffFwd.triangularView<Eigen::Lower>().solveInPlace(_coeffInv);
}

Eigen::MatrixXd HermiteTransformMatrix::compute(Eigen::Matrix2d const & transform, int order) const {
    if (order > _order) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            boost::str(
                boost::format("order passed to compute() (%d) is larger than construction order (%d)")
                % order % _order
            )
        );
    }
    int const size = computeSize(order);
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size, size);
    Binomial binomial_m(order);
    Binomial binomial_n(order);
    for (int jn=0, joff=0; jn <= order; joff += (++jn)) {
        for (int kn=jn, koff=joff; kn <= order; (koff += (++kn)) += (++kn)) {
            for (int jx=0,jy=jn; jx <= jn; ++jx,--jy) {
                for (int kx=0,ky=kn; kx <= kn; ++kx,--ky) {
                    double & element = result(koff+kx, joff+jx);
                    for (int m = 0; m <= order; ++m) {
                        int const order_minus_m = order - m;
                        binomial_m.reset(m, transform(0,0), transform(0,1));
                        for (int p = 0; p <= m; ++p) {
                            double const tmp1 = binomial_m[p] * _coeffFwd(kx, m);
                            for (int n = 0; n <= order_minus_m; ++n) {
                                binomial_n.reset(n, transform(1,0), transform(1,1));
                                double const tmp2 = _coeffFwd(ky, n) * tmp1;
                                for (int q = 0; q <= n; ++q) {
                                    element += tmp2 * _coeffInv(m+n-p-q, jx) * _coeffInv(p+q, jy)
                                        * binomial_n[q];
                                } // q
                            } // n
                        } // p
                    } // m
                } // kx,ky
            } // jx,jy
        } // kn
    } // jn
    return result;
}

}} // namespace lsst::shapelet
