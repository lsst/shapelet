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

#include "lsst/shapelet/ShapeletFunction.h"
#include "lsst/shapelet/GaussHermiteConvolution.h"
#include "lsst/afw/geom/Angle.h"
#include "ndarray/eigen.h"

namespace lsst { namespace shapelet {

//================= TripleProductIntegral ===================================================================

namespace {

class TripleProductIntegral {
public:

    ndarray::Array<double,3,3> const asArray() const { return _array; }

    static ndarray::Array<double,3,3> make1d(int const order1, int const order2, int const order3);

    TripleProductIntegral(int order1, int order2, int order3);

    ndarray::Vector<int,3> const & getOrders() const { return _orders; }

private:

    ndarray::Vector<int,3> _orders;
    ndarray::Array<double,3,3> _array;
};

TripleProductIntegral::TripleProductIntegral(int order1, int order2, int order3) :
    _orders(ndarray::makeVector(order1, order2, order3)),
    _array(
        ndarray::allocate(
            ndarray::makeVector(
                computeOffset(order1+1),
                computeOffset(order2+1),
                computeOffset(order3+1)
            )
        )
    )
{
    _array.deep() = 0.0;
    ndarray::Array<double,3,3> a1d = make1d(order1, order2, order3);
    ndarray::Array<double,3,3>::Index o, i, x, y, n;
    for (n[0] = o[0] = 0; n[0] <= _orders[0]; o[0] += ++n[0]) {
        for (n[2] = o[2] = 0; n[2] <= _orders[2]; o[2] += ++n[2]) {
            int n0n2 = n[0] + n[2];
            int max = (order2 < n0n2) ? _orders[1] : n0n2;
            for (o[1] = n[1] = bool(n0n2 % 2); n[1] <= max; (o[1] += ++n[1]) += ++n[1]) {
                for (i[0] = o[0], x[0] = 0, y[0] = n[0]; x[0] <= n[0]; ++x[0], --y[0], ++i[0]) {
                    for (i[2] = o[2], x[2] = 0, y[2] = n[2]; x[2] <= n[2]; ++x[2], --y[2], ++i[2]) {
                        for (i[1] = o[1], x[1] = 0, y[1] = n[1]; x[1] <= n[1]; ++x[1], --y[1], ++i[1]) {
                            _array[i] = a1d[x] * a1d[y];
                        }
                    }
                }
            }
        }
    }
}

ndarray::Array<double,3,3>
TripleProductIntegral::make1d(int order1, int order2, int order3) {
    ndarray::Array<double,3,3> array = ndarray::allocate(ndarray::makeVector(order1+1, order2+1, order3+1));
    array.deep() = 0.0;
    double first = array[ndarray::makeVector(0, 0, 0)] = BASIS_NORMALIZATION * M_SQRT1_2;
    if (order1 < 1) {
        if (order3 < 1) return array;
        // n1=0, n2=0, n3=even
        double previous = first;
        for (int n3 = 2; n3 <= order3; n3 += 2) {
            previous = array[ndarray::makeVector(0, 0, n3)] = -0.5*rationalSqrt(n3 - 1, n3)*previous;
        }
        // n1=0, n2>0, n3>0
        for (int n2 = 1; n2 <= order2; ++n2) {
            for (int n3 = n2; n3 <= order3; n3 += 2) {
                array[ndarray::makeVector(0, n2, n3)]
                    = rationalSqrt(n3, 2*n2) * array[ndarray::makeVector(0, n2-1, n3-1)];
            }
        }
        return array;
    }

    if (order3 < 1) {
        // n1=even, n2=0, n3=even
        double previous = first;
        for (int n1 = 2; n1 <= order1; n1 += 2) {
            previous = array[ndarray::makeVector(0, 0, n1)] = -0.5*rationalSqrt(n1 - 1, n1)*previous;
        }
        // n1=0, n2>0, n3>0
        for (int n2 = 1; n2 <= order2; ++n2) {
            for (int n1 = n2; n1 <= order1; n1 += 2) {
                array[ndarray::makeVector(n1, n2, 0)]
                    = rationalSqrt(n1, 2*n2) * array[ndarray::makeVector(n1-1, n2-1, 0)];
            }
        }
        return array;
    }

    // n1=any, n2=0, n3=(0,1)
    {
        double previous = first;
        int n1 = 0;
        while (true) {
            if (++n1 > order1) break;
            array[ndarray::makeVector(n1, 0, 1)] = 0.5*intSqrt(n1)*previous;
            if (++n1 > order1) break;
            previous = array[ndarray::makeVector(n1, 0, 0)] = -0.5*rationalSqrt(n1 - 1, n1)*previous;
        }
    }

    // n1=(0,1), n2=0, n3=any
    {
        double previous = first;
        int n3 = 0;
        while (true) {
            if (++n3 > order3) break;
            array[ndarray::makeVector(1, 0, n3)] = 0.5*intSqrt(n3) * previous;
            if (++n3 > order3) break;
            previous = array[ndarray::makeVector(0, 0, n3)] = -0.5*rationalSqrt(n3 - 1, n3)*previous;
        }
    }

    // n1>1, n2=0, n3>1
    for (int n1 = 2; n1 <= order1; ++n1) {
        double f1 = -0.5*rationalSqrt(n1 - 1, n1);
        for (int n3 = (n1 % 2) ? 3:2; n3 <= order3; n3 += 2) {
            array[ndarray::makeVector(n1, 0, n3)]
                = f1 * array[ndarray::makeVector(n1-2, 0, n3)]
                + 0.5*rationalSqrt(n3, n1) * array[ndarray::makeVector(n1-1, 0, n3-1)];
        }
    }

    for (int n2 = 1; n2 <= order2; ++n2) {
        // n1>=n2, n3=0
        for (int n1 = n2; n1 <= order1; n1 += 2) {
            array[ndarray::makeVector(n1, n2, 0)]
                = rationalSqrt(n1, 2*n2) * array[ndarray::makeVector(n1-1, n2-1, 0)];
        }
        // n1=0, n3>=n2
        for (int n3 = n2; n3 <= order3; n3 += 2) {
            array[ndarray::makeVector(0, n2, n3)]
                = rationalSqrt(n3, 2*n2) * array[ndarray::makeVector(0, n2-1, n3-1)];
        }
        // n1>0, n3>0
        for (int n1 = 1; n1 <= order1; ++n1) {
            for (int n3 = ((n1+n2) % 2) ? 1:2; n3 <= order3; n3 += 2) {
                array[ndarray::makeVector(n1, n2, n3)]
                    = rationalSqrt(n1, 2*n2) * array[ndarray::makeVector(n1-1, n2-1, n3)]
                    + rationalSqrt(n3, 2*n2) * array[ndarray::makeVector(n1, n2-1, n3-1)];
            }
        }
    }

    return array;
}

} // anonymous

//================= Impl Base Class =========================================================================

class GaussHermiteConvolution::Impl {
public:

    Impl(int colOrder, ShapeletFunction const & psf);

    virtual ndarray::Array<double const,2,2> evaluate(afw::geom::ellipses::Ellipse & ellipse) const = 0;

    int getColOrder() const { return _colOrder; }

    int getRowOrder() const { return _rowOrder; }

    virtual ~Impl() {}

protected:
    int _rowOrder;
    int _colOrder;
    ShapeletFunction _psf;
    ndarray::Array<double,2,2> _result;
};

GaussHermiteConvolution::Impl::Impl(
    int colOrder, ShapeletFunction const & psf
) :
    _rowOrder(colOrder + psf.getOrder()), _colOrder(colOrder), _psf(psf),
    _result(ndarray::allocate(computeSize(_rowOrder), computeSize(_colOrder)))
{
    _psf.changeBasisType(HERMITE);
}

//================= ImplN: Arbitrary Order ==================================================================

namespace {

class ImplN : public GaussHermiteConvolution::Impl {
public:

    ImplN(int colOrder, ShapeletFunction const & psf);

    virtual ndarray::Array<double const,2,2> evaluate(afw::geom::ellipses::Ellipse & ellipse) const;

private:
    TripleProductIntegral _tpi;
    HermiteTransformMatrix _htm;
};

ImplN::ImplN(
    int colOrder, ShapeletFunction const & psf
) :
    GaussHermiteConvolution::Impl(colOrder, psf),
    _tpi(psf.getOrder(), _rowOrder, _colOrder),
    _htm(_rowOrder)
{}

ndarray::Array<double const,2,2> ImplN::evaluate(
    afw::geom::ellipses::Ellipse & ellipse
) const {
    auto result = ndarray::asEigenMatrix(_result);
    auto psf_coeff = ndarray::asEigenMatrix(_psf.getCoefficients());

    Eigen::Matrix2d psfT = _psf.getEllipse().getCore().getGridTransform().invert().getMatrix();
    Eigen::Matrix2d modelT = ellipse.getCore().getGridTransform().invert().getMatrix();
    ellipse.convolve(_psf.getEllipse()).inPlace();
    Eigen::Matrix2d convolvedTI = ellipse.getCore().getGridTransform().getMatrix() * std::sqrt(2.0);
    Eigen::Matrix2d psfArg = (convolvedTI * psfT).transpose();
    Eigen::Matrix2d modelArg = (convolvedTI * modelT).transpose();

    int const psfOrder = _psf.getOrder();

    Eigen::MatrixXd psfMat = _htm.compute(psfArg, psfOrder);
    Eigen::MatrixXd modelMat = _htm.compute(modelArg, _colOrder);

    // [kq]_m = \sum_m i^{n+m} [psfMat]_{m,n} [psf]_n
    // kq is zero unless {n+m} is even
    Eigen::VectorXd kq = Eigen::VectorXd::Zero(psfMat.size());
    for (int m = 0, mo = 0; m <= psfOrder; mo += ++m) {
        for (int n = m, no = mo; n <= psfOrder; no += ++n, no += ++n) {
            if ((n + m) % 4) { // (n + m) % 4 is always 0 or 2
                kq.segment(mo, m+1) -= psfMat.block(no, mo, n+1, m+1).adjoint() * psf_coeff.segment(no, n+1);
            } else {
                kq.segment(mo, m+1) += psfMat.block(no, mo, n+1, m+1).adjoint() * psf_coeff.segment(no, n+1);
            }
        }
    }
    kq *= 4.0 * afw::geom::PI;

    // [kqb]_{m,n} = \sum_l i^{m-n-l} [kq]_l [tpi]_{l,m,n}
    Eigen::MatrixXd kqb = Eigen::MatrixXd::Zero(result.rows(), result.cols());
    {
        ndarray::Array<double,3,3> b = _tpi.asArray();
        ndarray::Vector<int,3> n, o, x;
        x[0] = 0;
        for (n[1] = o[1] = 0; n[1] <= _rowOrder; o[1] += ++n[1]) {
            for (n[2] = o[2] = 0; n[2] <= _colOrder; o[2] += ++n[2]) {
                for (n[0] = o[0] = bool((n[1]+n[2])%2); n[0] <= psfOrder; o[0] += ++n[0], o[0] += ++n[0]) {
                    if (n[0] + n[2] < n[1]) continue; // b is triangular
                    double factor = ((n[2] - n[0] - n[1]) % 4) ? -1.0 : 1.0;
                    for (x[1] = 0; x[1] <= n[1]; ++x[1]) {
                        for (x[2] = 0; x[2] <= n[2]; ++x[2]) {
                            // TODO: optimize this inner loop - it's sparse
                            double & out = kqb(o[1] + x[1], o[2] + x[2]);
                            for (x[0] = 0; x[0] <= n[0]; ++x[0]) {
                                out += factor * kq[o[0] + x[0]] * b[o + x];
                            }
                        }
                    }
                }
            }
        }
    }

    result.setZero();
    for (int m = 0, mo = 0; m <= _rowOrder; mo += ++m) {
        for (int n = 0, no = 0; n <= _colOrder; no += ++n) {
            int jo = bool(n % 2);
            for (int j = jo; j <= n; jo += ++j, jo += ++j) {
                if ((n - j) % 4) { // (n - j) % 4 is always 0 or 2
                    result.block(mo, no, m+1, n+1)
                        -= kqb.block(mo, jo, m+1, j+1) * modelMat.block(no, jo, n+1, j+1).adjoint();
                } else {
                    result.block(mo, no, m+1, n+1)
                        += kqb.block(mo, jo, m+1, j+1) * modelMat.block(no, jo, n+1, j+1).adjoint();
                }
            }
        }
    }

    return _result;
}

} // anonymous

//================= Impl0: Zeroth Order =====================================================================

namespace {

class Impl0 : public GaussHermiteConvolution::Impl {
public:

    typedef double Scalar;

    explicit Impl0(ShapeletFunction const & psf);

    virtual ndarray::Array<double const,2,2> evaluate(afw::geom::ellipses::Ellipse & ellipse) const;

private:
    ndarray::Array<Scalar,1,0> _r;
    ndarray::Array<Scalar,4,4> _p;
};

Impl0::Impl0(
    ShapeletFunction const & psf
) :
    GaussHermiteConvolution::Impl(0, psf),
    _r(_result[ndarray::view()(0)]), // _result has one column, so it's more convenient if we have a 1-d view
    _p(ndarray::allocate(ndarray::makeVector(_rowOrder+1, _rowOrder+1, _rowOrder+1, _rowOrder+1)))
{
    _p.deep() = 0.0;
}

ndarray::Array<double const,2,2> Impl0::evaluate(
    afw::geom::ellipses::Ellipse & ellipse
) const {
    ndarray::Array<double const,1,1> psfCoeff(_psf.getCoefficients());

    Eigen::Matrix2d psfT = _psf.getEllipse().getCore().getGridTransform().invert().getMatrix();
    ellipse.convolve(_psf.getEllipse()).inPlace();
    Eigen::Matrix2d convolvedTI = ellipse.getCore().getGridTransform().getMatrix();
    Eigen::Matrix2d psfArg = (convolvedTI * psfT).transpose();

    int const psfOrder = _psf.getOrder();

    // Just for readability: extract matrix elements and strides
    double const sxx = psfArg(0,0);
    double const sxy = psfArg(0,1);
    double const syx = psfArg(1,0);
    double const syy = psfArg(1,1);

    int const idx_00_01 = _p.getStride<3>();
    int const idx_00_10 = _p.getStride<2>();
    int const idx_01_00 = _p.getStride<1>();
    int const idx_10_00 = _p.getStride<0>();

    // at each step in many of the loops, we increment n_x and decrement n_y (or likewise with m)
    int const incr_n = idx_10_00 - idx_01_00;
    int const incr_m = idx_00_10 - idx_00_01;

    //---- Some terminology to help with the ugly math and abbreviated comments below

    // The offsetX variables below are the offsets of the elements on the rhs of
    // the recurrence relations, relative to the element we're setting on the
    // lhs.  Hence a recurrence setting line will often look like this:
    //
    //    *pc = t1 * pc[offset1] + t2 * pc[offset2];
    //
    // Here, "pc" is a pointer to the current item, so we set it simply by
    // dereferencing it, and we refer to the previous elements using the offsetX
    // variables.  Number suffixes indicate the term in the recurrence, while
    // letter suffices indicate different related recurrence relations.

    double * p = _p.getData();
    *p = 1.0;

    // do n == 0, m >= 2, m even
    if (psfOrder > 1) {
        double const t1a = sxx*sxx + sxy*sxy - 1;
        double const t1b = syy*syy + syx*syx - 1;
        double const t2 = sxx*syx + sxy*syy;
        int const offset1a = -2*idx_00_10;
        int const offset1b = -2*idx_00_01;
        int const offset2 = -idx_00_10 - idx_00_01;
        for (int m = 2; m <= psfOrder; m += 2) {
            int m_x = 0, m_y = m;
            double * pc = p + idx_00_01*m_y;
            // for m_x == 0, use 'b' recurrence, dropping terms that vanish for m_x == 0:
            // p[n_x,n_y,m_x,m_y] = \sqrt{\frac{m_y-1}{m_y}} (S_{yy}^2 + S_{yx}^2 - 1) p[n_x,n_y,m_x,m_y-2]
            for (; m_x == 0; ++m_x, --m_y, pc += incr_m) {
                *pc = t1b*rationalSqrt(m_y-1,m_y)*pc[offset1b];
            }
            // for m_x == 1, use 'a' recurrence, dropping first term
            for (; m_x == 1; ++m_x, --m_y, pc += incr_m) {
                *pc = t2*rationalSqrt(m_y, m_x)*pc[offset2];
            }
            // for m_x > 0, m_y > 0, use full 'a' recurrence:
            // p[n_x,n_y,m_x,m_y] = \sqrt{\frac{m_x-1}{m_x}} (S_{xx}^2 + S_{xy}^2 - 1) p[n_x,n_y,m_x-2,m_y]
            //                    + \sqrt{m_y/m_x} (S_{xx}S_{yx} + S_{xy}S_{yy}) p[n_x,n_y,m_x-1,m_y-1]
            for (; m_y > 0; ++m_x, --m_y, pc += incr_m) {
                *pc = t1a*rationalSqrt(m_x-1,m_x)*pc[offset1a]
                    + t2*rationalSqrt(m_y, m_x)*pc[offset2];
            }
            // for m_y == 0, use 'a' recurrence, drop 2nd term
            for (; m_y == 0; ++m_x, --m_y, pc += incr_m) {
                *pc = t1a*rationalSqrt(m_x-1,m_x)*pc[offset1a];
            }
        }
    }

    // do n > 0, m >= n, m + n even
    if (psfOrder > 0) {
        int const offset1a = -idx_10_00 - idx_00_10;
        int const offset2a = -idx_10_00 - idx_00_01;
        int const offset1b = -idx_01_00 - idx_00_01;
        int const offset2b = -idx_01_00 - idx_00_10;
        for (int n = 1; n <= psfOrder; ++n) {
            int n_x = 0, n_y = n;
            double * pnc = p + idx_01_00*n_y;
            // for n_x == 0, use 'b' recurrence
            // p[n_x,n_y,m_x,m_y] = \sqrt{m_x/n_x} S_{xx} p[n_x-1,n_y,m_x-1,m_y]
            //                    + \sqrt{m_y,n_x} S_{yx} p[n_x-1,n_y,m_x,m_y-1]
            for (; n_x == 0; ++n_x, --n_y, pnc += incr_n) {
                for (int m = n; m <= psfOrder; m += 2) {
                    int m_x = 0, m_y = m;
                    double * pc = pnc + idx_00_01*m_y;
                    // for m_x == 0, drop 2nd term
                    for (; m_x == 0; ++m_x, --m_y, pc += incr_m) {
                        *pc = rationalSqrt(m_y, n_y)*syy*pc[offset1b];
                    }
                    // for m_x > 0, m_y > 0, use both terms
                    for (; m_x < m; ++m_x, --m_y, pc += incr_m) {
                        *pc = rationalSqrt(m_y, n_y)*syy*pc[offset1b]
                            + rationalSqrt(m_x, n_y)*sxy*pc[offset2b];
                    }
                    // for m_y == 0, drop 1st term
                    for (; m_y == 0; ++m_x, --m_y, pc += incr_m) {
                        *pc = rationalSqrt(m_x, n_y)*sxy*pc[offset2b];
                    }
                }
            }
            // for n_x > 0, use 'a' recurrence:
            // p[n_x,n_y,m_x,m_y] = \sqrt{m_y/n_y} S_{yy} p[n_x,n_y-1,m_x,m_y-1]
            //                    + \sqrt{m_x,n_y} S_{xy} p[n_x,n_y-1,m_x-1,m_y]
            for (; n_x <= n; ++n_x, --n_y, pnc += incr_n) {
                for (int m = n; m <= psfOrder; m += 2) {
                    int m_x = 0, m_y = m;
                    double * pc = pnc + idx_00_01*m_y;
                    // for m_x == 0, drop 1st term
                    for (; m_x == 0; ++m_x, --m_y, pc += incr_m) {
                        *pc = rationalSqrt(m_y, n_x)*syx*pc[offset2a];
                    }
                    // for m_x > 0, m_y > 0, use both terms
                    for (; m_x < m; ++m_x, --m_y, pc += incr_m) {
                        *pc = rationalSqrt(m_x, n_x)*sxx*pc[offset1a]
                            + rationalSqrt(m_y, n_x)*syx*pc[offset2a];
                    }
                    // for m_y == 0, drop 2nd term
                    for (; m_y == 0; ++m_x, --m_y, pc += incr_m) {
                        *pc = rationalSqrt(m_x, n_x)*sxx*pc[offset1a];
                    }
                }
            }
        }
    }

    // Form inner product between p and psfCoeffs, including 2\sqrt{\pi} i^{m-n} factor
    double const f1 = 2.0 * std::sqrt(M_PI);
    for (int n = 0, n_o = 0; n <= psfOrder; n_o += ++n) {
        double const * pnc = p + idx_01_00*n;
        ndarray::Array<double,1,0>::Iterator rIter = _r.begin() + n_o;
        for (int n_x = 0; n_x <= n; ++n_x, ++rIter, pnc += incr_n) {
            *rIter = 0;
            for (int m = n, m_o = computeSize(n-1); m <= psfOrder; m_o += ++m, m_o += ++m) {
                double f2 = ((m-n) % 4) ? -f1 : f1;
                double const * pc = pnc + idx_00_01*m;
                ndarray::Array<double const,1,1>::Iterator psfIter = psfCoeff.begin() + m_o;
                for (int m_x = 0; m_x <= m; ++m_x, ++psfIter, pc += incr_m) {
                    *rIter += f2 * (*pc) * (*psfIter);
                }
            }
        }
    }

    return _result;
}

} // anonymous

//================= public GaussHermiteConvolution ==========================================================

int GaussHermiteConvolution::getRowOrder() const { return _impl->getRowOrder(); }

int GaussHermiteConvolution::getColOrder() const { return _impl->getColOrder(); }

ndarray::Array<double const,2,2>
GaussHermiteConvolution::evaluate(
    afw::geom::ellipses::Ellipse & ellipse
) const {
    return _impl->evaluate(ellipse);
}

GaussHermiteConvolution::GaussHermiteConvolution(
    int colOrder,
    ShapeletFunction const & psf
) : _impl() {
    if (colOrder == 0) {
        _impl.reset(new Impl0(psf));
    } else {
        _impl.reset(new ImplN(colOrder, psf));
    }
}

GaussHermiteConvolution::~GaussHermiteConvolution() {}

}} // namespace lsst::shapelet
