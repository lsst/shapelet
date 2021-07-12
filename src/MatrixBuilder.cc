// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#include <cmath>
#include <memory>
#include "ndarray/eigen.h"

#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/shapelet/GaussHermiteConvolution.h"

namespace lsst { namespace shapelet {

//===========================================================================================================
//================== Implementation Base Classes ============================================================
//===========================================================================================================

/*
 * The implementation classes held by MatrixBuilder and MatrixBuilderFactory are actually just abstract
 * base classes - the concrete derived classes (in anonymous namespaces, below) correspond to different
 * constructor arguments.
 *
 * For each Builder/Factory implementation pair, we declare the Factory as an inner class of the Builder
 * class, just to avoid template and forward-declaration boilerplate.
 *
 * The public MatrixBuilderFactory class holds a MatrixBuilderFactory::Impl, which it uses to create a
 * MatrixBuilder::Impl, which is passed to the MatrixBuilder private constructor.  When one of the
 * MatrixBuilder public constructors is called directly, it simply forwards those arguments to
 * a MatrixBuilderFactory constructor, then uses the MatrixBuilderFactory to create the
 * MatrixBuilder::Impl it will hold.
 *
 * Which Impl derived classes are used is set by the arguments to the MatrixBuilderFactory or
 * MatrixBuilder constructor.
 */

template <typename T>
class MatrixBuilder<T>::Impl {
public:

    virtual int getDataSize() const = 0;

    virtual int getBasisSize() const = 0;

    virtual void buildMatrix(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) = 0;

    virtual ~Impl() {}

};

template <typename T>
class MatrixBuilderFactory<T>::Impl {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;
    typedef typename MatrixBuilder<T>::Impl BuilderImpl;

    virtual int getDataSize() const = 0;

    virtual int getBasisSize() const = 0;

    virtual int computeWorkspace() const = 0;

    std::shared_ptr<BuilderImpl> makeBuilderImpl() const {
        Workspace workspace(computeWorkspace());
        return makeBuilderImpl(workspace);
    }

    virtual std::shared_ptr<BuilderImpl> makeBuilderImpl(Workspace & workspace) const = 0;

    virtual ~Impl() {}

};

//===========================================================================================================
//================== Single-Component Implementation Base Classes ===========================================
//===========================================================================================================

/*
 * This implementation pair is an intermediate base class for cases where we have only one shapelet basis.
 * It holds the original coordinate arrays and a set of temporary arrays that are filled with the values
 * of the coordinate arrays in the coordinate system where the ellipse we're evaluating on is the unit
 * circle.  In other words, all the derived classes just have to worry about evaluating circular basis
 * functions, and the base class handles turning those into ellipses by giving them transformed inputs.
 */

namespace {

template <typename T>
class SimpleImpl : public MatrixBuilder<T>::Impl {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    class Factory : public MatrixBuilderFactory<T>::Impl {
    public:

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y
        ) : _x(x), _y(y) {}

        virtual int getDataSize() const { return _x.template getSize<0>(); }

        virtual int computeWorkspace() const { return 2*_x.template getSize<0>(); }

        ndarray::Array<T const,1,1> getX() const { return _x; }

        ndarray::Array<T const,1,1> getY() const { return _y; }

    private:
        ndarray::Array<T const,1,1> _x;
        ndarray::Array<T const,1,1> _y;
    };

    SimpleImpl(Factory const & factory, Workspace * workspace) :
        _x(factory.getX()), _y(factory.getY()),
        _xt(workspace->makeVector(factory.getDataSize())),
        _yt(workspace->makeVector(factory.getDataSize())),
        _detFactor(1.0),
        _manager(workspace->getManager())
    {}

    virtual int getDataSize() const { return _x.template getSize<0>(); }

    void readEllipse(afw::geom::ellipses::Ellipse const & ellipse) {
        geom::AffineTransform transform = ellipse.getGridTransform();
        _xt = ndarray::asEigenArray(_x) * transform[geom::AffineTransform::XX]
            + ndarray::asEigenArray(_y) * transform[geom::AffineTransform::XY]
            + transform[geom::AffineTransform::X];
        _yt = ndarray::asEigenArray(_x) * transform[geom::AffineTransform::YX]
            + ndarray::asEigenArray(_y) * transform[geom::AffineTransform::YY]
            + transform[geom::AffineTransform::Y];
        _detFactor = transform.getLinear().computeDeterminant();
    }

private:
    ndarray::Array<T const,1,1> _x;
    ndarray::Array<T const,1,1> _y;
protected:
    Eigen::Map< Eigen::Array<T,Eigen::Dynamic,1> > _xt;
    Eigen::Map< Eigen::Array<T,Eigen::Dynamic,1> > _yt;
    T _detFactor;
private:
    ndarray::Manager::Ptr _manager;
};

} // anonymous

//===========================================================================================================
//================== Non-Convolved, Non-Remapped Shapelet Implementation ====================================
//===========================================================================================================

/*
 * This implementation pair is for a normal Gauss-Hermite basis, and it does the work of evaluating that
 * basis using 1-d recurrence relations for Hermite polynomials on x and y workspaces (we do all points at
 * once to take advantage of Eigen's vectorized math), and then zips these together and multiplies them
 * by the Gaussian term (also evaluated using vectorized Eigen math) to form the 2-d shapelet functions.
 *
 * For more information on the definition of shapelet functions and the ordering of the packing, see
 * the description of the HERMITE enum value in constants.h, or the PackedIndex class in
 * GaussHermiteEvaluator.
 *
 * The matrix produced by this implementation is used by all the other SimpleImpl derived classes below
 * (in fact, they all inherit directly or indirectly from ShapeletImpl), and in all of those cases it's
 * the left-hand side of a matrix equation.  To avoid confusion later, we'll refer to the order of the
 * shapelet basis produced by ShapeletImpl as the "lhs order".
 */

namespace {

template <typename T>
class ShapeletImpl : public SimpleImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    class Factory : public SimpleImpl<T>::Factory {
    public:

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y,
            int lhsOrder
        ) : SimpleImpl<T>::Factory(x, y),
            _lhsOrder(lhsOrder)
        {}

        virtual int getBasisSize() const { return computeSize(_lhsOrder); }

        virtual int computeWorkspace() const {
            return this->getDataSize()*(1 + 2*(_lhsOrder + 1))
                + SimpleImpl<T>::Factory::computeWorkspace();
        }

        int getLhsOrder() const { return _lhsOrder; }

        virtual std::shared_ptr<typename MatrixBuilder<T>::Impl> makeBuilderImpl(Workspace & workspace) const {
            return std::make_shared<ShapeletImpl>(*this, &workspace);
        }

    private:
        int _lhsOrder;
    };

    ShapeletImpl(Factory const & factory, Workspace * workspace) :
        SimpleImpl<T>(factory, workspace),
        _lhsOrder(factory.getLhsOrder()),
        _gaussian(workspace->makeVector(factory.getDataSize())),
        _xHermite(workspace->makeMatrix(factory.getDataSize(), factory.getLhsOrder() + 1)),
        _yHermite(workspace->makeMatrix(factory.getDataSize(), factory.getLhsOrder() + 1))
    {}

    virtual int getBasisSize() const { return computeSize(_lhsOrder); }

    virtual void buildMatrix(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        buildMatrix(ndarray::asEigenArray(output), ellipse);
    }

    template <typename EigenArrayT>
    void buildMatrix(
        EigenArrayT output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        this->readEllipse(ellipse);
        fillGaussian();
        fillHermite1d(_xHermite, this->_xt);
        fillHermite1d(_yHermite, this->_yt);
        for (PackedIndex i; i.getOrder() <= _lhsOrder; ++i) {
            output.col(i.getIndex())
                += this->_detFactor * _gaussian * _xHermite.col(i.getX()) * _yHermite.col(i.getY());
        }
    }

    void fillGaussian() {
        _gaussian = (-0.5*(this->_xt.square() + this->_yt.square())).exp();
    }

    template <typename CoordArray>
    void fillHermite1d(
        typename Workspace::Matrix & output,
        CoordArray const & coord
    ) {
        if (output.cols() > 0) {
            output.col(0).setConstant(BASIS_NORMALIZATION);
        }
        if (output.cols() > 1) {
            output.col(1) = intSqrt(2) * coord * output.col(0);
        }
        for (int j = 2; j <= _lhsOrder; ++j) {
            output.col(j) = rationalSqrt(2, j) * coord * output.col(j-1)
                - rationalSqrt(j - 1, j) * output.col(j-2);
        }
    }

protected:
    int _lhsOrder;
    typename Workspace::Vector _gaussian;
    typename Workspace::Matrix _xHermite;
    typename Workspace::Matrix _yHermite;
};

} // anonymous


//===========================================================================================================
//================== Convolved, Non-Remapped Shapelet Implementation ========================================
//===========================================================================================================

/*
 * This implementation pair handles a Gauss-Hermite basis convolved with a ShapeletFunction.  To convolve
 * shapelets, we need to do three things:
 *
 *  - We need to convolve the basis ellipse with the PSF ellipse, and use that to evaluate the actual
 *    shapelet basis (by "convolving ellipses", I actually mean we imagine the input ellipses as
 *    corresponding to a pair of elliptical Gaussians, and convolve those, and get the ellipse of the
 *    result).
 *
 *  - We need to evaluate a regular shapelet basis, using the "convolved" ellipse, at a higher order
 *    than the convolved basis we want.  In the code, we call the higher-order basis matrix the "lhs"
 *    matrix, and this higher order the "lhs order".
 *
 *  - We multiply the lhs basis matrix by another matrix that maps the lhs basis to the
 *    lower-order rhs basis that corresponds to the requested order of the total.  Note that the
 *    lhs matrix has dimensions (dataSize, computeSize(lhsOrder)), while the rhs convolution matrix
 *    has dimensions (computeSize(lhsOrder), computeSize(rhsOrder)), so their product has dimensions
 *    (dataSize, computeSize(rhsOrder)), which is of course the expected size of the output.
 *
 * So, the procedure for evaluating a convolved shapelet basis looks like this:
 *  1) call GaussianHermiteConvolution::evaluate() to modify the ellipse and get back the convolution matrix
 *  2) evaluate a standard shapelet basis using the new ellipse, and store it in a temporary (this is the
 *     workspace matrix "_lhs").  This is delegated to the ShapeletImpl base class.
 *  3) compute the matrix product.
 */

namespace {

template <typename T>
class ConvolvedShapeletImpl : public ShapeletImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    class Factory : public ShapeletImpl<T>::Factory {
    public:

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y,
            int rhsOrder,
            ShapeletFunction const & psf
        ) :
            ShapeletImpl<T>::Factory(
                x, y, GaussHermiteConvolution::computeRowOrder(rhsOrder, psf.getOrder())
            ),
            _rhsOrder(rhsOrder),
            _psf(psf)
        {}

        virtual int getBasisSize() const { return computeSize(_rhsOrder); }

        virtual int computeWorkspace() const {
            return ShapeletImpl<T>::Factory::computeWorkspace()
                +  this->getDataSize() * computeSize(this->getLhsOrder());
        }

        int getRhsOrder() const { return _rhsOrder; }

        ShapeletFunction const & getPsf() const { return _psf; }

        virtual std::shared_ptr<typename MatrixBuilder<T>::Impl> makeBuilderImpl(Workspace & workspace) const {
            return std::make_shared<ConvolvedShapeletImpl>(*this, &workspace);
        }

    protected:
        int _rhsOrder;
        ShapeletFunction _psf;
    };

    ConvolvedShapeletImpl(Factory const & factory, Workspace * workspace) :
        ShapeletImpl<T>(factory, workspace),
        _ellipse(afw::geom::ellipses::Quadrupole()),
        _convolution(factory.getRhsOrder(), factory.getPsf()),
        _lhs(workspace->makeMatrix(factory.getDataSize(), computeSize(_convolution.getRowOrder())))
    {}

    virtual int getBasisSize() const { return computeSize(_convolution.getColOrder()); }

    virtual void buildMatrix(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        if (!(ellipse.getCore().getDeterminantRadius() >= 0.0)) {
            throw LSST_EXCEPT(
                pex::exceptions::UnderflowError,
                "Underflow error in ellipse scaling/convolution"
            );
        }
        _ellipse = ellipse;
        ndarray::Array<double const,2,2> rhs = computeTerms();
        ndarray::asEigenMatrix(output) += _lhs.matrix() * ndarray::asEigenMatrix(rhs).cast<T>();
    }

    ndarray::Array<double const,2,2> computeTerms() {
        ndarray::Array<double const,2,2> rhs = _convolution.evaluate(_ellipse);
        _lhs.setZero();
        ShapeletImpl<T>::buildMatrix(_lhs, _ellipse);
        return rhs;
    }

protected:
    mutable afw::geom::ellipses::Ellipse _ellipse;
    GaussHermiteConvolution _convolution;
    typename Workspace::Matrix _lhs;
};

} // anonymous


//===========================================================================================================
//================== Non-Convolved, Remapped Shapelet Implementation ========================================
//===========================================================================================================

/*
 * This implementation pair handles a standard Gauss-Hermite basis that is remapped to a different
 * basis via a "remap" matrix.  That is, each of the elements of the remapped matrix is a linear
 * combination of the elements in the original Gauss-Hermite basis.
 *
 * The procedure for evaluating such a basis is pretty straightforward: we evaluate the standard
 * Gauss-Hermite basis matrix using the ShapeletImpl derived class, store that in a temporary (again
 * called _lhs), and then multiply that by the remap matrix to get the full remapped basis matrix.
 */

namespace {

template <typename T>
class RemappedShapeletImpl : public ShapeletImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    class Factory : public ShapeletImpl<T>::Factory {
    public:

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y,
            int lhsOrder,
            double radius,
            ndarray::Array<double const,2,2> const & remapMatrix
        ) : ShapeletImpl<T>::Factory(x, y, lhsOrder),
            _radius(radius),
            _remapMatrix(remapMatrix)
        {}

        virtual int getBasisSize() const { return _remapMatrix.getSize<1>(); }

        virtual int computeWorkspace() const {
            return ShapeletImpl<T>::Factory::computeWorkspace()
                +  this->getDataSize() * computeSize(this->getLhsOrder());
        }

        double getRadius() const { return _radius; }

        ndarray::Array<double const,2,2> getRemapMatrix() const { return _remapMatrix; }

        virtual std::shared_ptr<typename MatrixBuilder<T>::Impl> makeBuilderImpl(Workspace & workspace) const {
            return std::make_shared<RemappedShapeletImpl>(*this, &workspace);
        }

    protected:
        double _radius;
        ndarray::Array<double const,2,2> _remapMatrix;
    };

    RemappedShapeletImpl(
        Factory const & factory,
        Workspace * workspace
    ) : ShapeletImpl<T>(factory, workspace),
        _ellipse(afw::geom::ellipses::Quadrupole()),
        _radius(factory.getRadius()),
        // transpose the remap matrix to preserve memory order when we copy it; will untranspose later
        _remapMatrix(ndarray::asEigenMatrix(factory.getRemapMatrix()).template cast<T>().transpose()),
        _lhs(workspace->makeMatrix(factory.getDataSize(), computeSize(factory.getLhsOrder())))
    {}

    virtual int getBasisSize() const { return _remapMatrix.rows(); }

    virtual void buildMatrix(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        _ellipse = ellipse;
        _ellipse.scale(_radius);
        _lhs.setZero();
        ShapeletImpl<T>::buildMatrix(_lhs, _ellipse);
        // undo the transpose in the constructor
        ndarray::asEigenMatrix(output) += _lhs.matrix() * _remapMatrix.transpose();
    }

protected:
    mutable afw::geom::ellipses::Ellipse _ellipse;
    double _radius;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> _remapMatrix;
    typename Workspace::Matrix _lhs;
};

} // anonymous


//===========================================================================================================
//================== Convolved, Remapped Shapelet Implementation ========================================
//===========================================================================================================

/*
 * This implementation pair is a combination of the convolved and remapped implementations above (we
 * inherit from the convolved implementation, and essentially reimplement the remapping, as there's
 * not much we can reuse in practice).  In this case, we have three terms in the matrix product:
 *
 *   - the "lhs" standard shapelet basis, evaluated at a higher order and a modified ellipse, as
 *     described above in the convolved implementation.
 *
 *   - the convolution matrix
 *
 *   - the remap matrix
 *
 * Because the convolution matrix and remap matrix have smaller dimensions than the lhs matrix, it's
 * generally more efficient to multiply the latter two terms first (we call the result the "rhs" here),
 * and then multiply the lhs matrix with that.
 */

namespace {

template <typename T>
class RemappedConvolvedShapeletImpl : public ConvolvedShapeletImpl<T> {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    class Factory : public ConvolvedShapeletImpl<T>::Factory {
    public:

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y,
            int rhsOrder,
            double radius,
            ndarray::Array<double const,2,2> const & remapMatrix,
            ShapeletFunction const & psf
        ) : ConvolvedShapeletImpl<T>::Factory(x, y, rhsOrder, psf),
            _radius(radius),
            _remapMatrix(remapMatrix)
        {}

        virtual int getBasisSize() const { return _remapMatrix.getSize<1>(); }

        virtual int computeWorkspace() const {
            return ConvolvedShapeletImpl<T>::Factory::computeWorkspace()
                + computeSize(this->getLhsOrder()) * getBasisSize();
        }

        double getRadius() const { return _radius; }

        ndarray::Array<double const,2,2> getRemapMatrix() const { return _remapMatrix; }

        virtual std::shared_ptr<typename MatrixBuilder<T>::Impl> makeBuilderImpl(Workspace & workspace) const {
            return std::make_shared<RemappedConvolvedShapeletImpl>(*this, &workspace);
        }

    protected:
        double _radius;
        ndarray::Array<double const,2,2> _remapMatrix;
    };

    RemappedConvolvedShapeletImpl(Factory const & factory, Workspace * workspace) :
        ConvolvedShapeletImpl<T>(factory, workspace),
        _radius(factory.getRadius()),
        // transpose the remap matrix to preserve memory order when we copy it; will untranspose later
        _remapMatrix(ndarray::asEigenMatrix(factory.getRemapMatrix()).template cast<T>().transpose()),
        _rhs(workspace->makeMatrix(computeSize(factory.getLhsOrder()), factory.getBasisSize()))
    {}

    virtual int getBasisSize() const { return _remapMatrix.rows(); }

    virtual void buildMatrix(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        this->_ellipse = ellipse;
        this->_ellipse.scale(_radius);
        ndarray::Array<double const,2,2> convolutionMatrix = this->computeTerms();
        _rhs.matrix() = ndarray::asEigenMatrix(convolutionMatrix).cast<T>() * _remapMatrix.transpose();  // untranspose
        ndarray::asEigenMatrix(output) += this->_lhs.matrix() * _rhs.matrix();
    }

protected:
    double _radius;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> _remapMatrix;
    typename Workspace::Matrix _rhs;
};

} // anonymous


//===========================================================================================================
//================== Multi-Component Implementations ========================================================
//===========================================================================================================

/*
 * This implementation pair is used for all cases where we have multiple shapelet components - those with
 * a MultiShapeletBasis and/or MultiShapeletFunction psf on input.
 *
 * In this case, the factory is simply a vector of single-component factories, and the builder is just a
 * vector of builders.  We always use the remapped single-component implementations, because a
 * MultiShapeletBasis has to have a non-identity remap matrix if there's more than one component (and if
 * there is just one component, the public class constructors are smart enough to just use a single-component
 * implementation).
 *
 * The only tricky bit here is the workspace management: instead of letting each component have a piece
 * of a large workspace, the components actually share the same workspace, so the workspace needed by
 * the compound implementation is actually the maximum needed by any of its components, not the sum of the
 * workspace needed by its components.
 */

namespace {

template <typename T>
class CompoundImpl : public MatrixBuilder<T>::Impl {
public:

    typedef MatrixBuilderWorkspace<T> Workspace;

    typedef typename MatrixBuilder<T>::Impl Component;
    typedef std::vector<std::shared_ptr<Component>> Vector;
    typedef typename Vector::const_iterator Iterator;

    typedef typename MatrixBuilderFactory<T>::Impl FactoryComponent;
    typedef std::vector<std::shared_ptr<FactoryComponent>> FactoryVector;
    typedef typename FactoryVector::const_iterator FactoryIterator;

    class Factory : public MatrixBuilderFactory<T>::Impl {
    public:

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y,
            MultiShapeletBasis const & basis
        ) {
            _components.reserve(basis.getComponentCount());
            for (MultiShapeletBasis::Iterator i = basis.begin(); i != basis.end(); ++i) {
                _components.push_back(
                    std::make_shared< typename RemappedShapeletImpl<T>::Factory >(
                        x, y, i->getOrder(), i->getRadius(), i->getMatrix()
                    )
                );
            }
        }

        Factory(
            ndarray::Array<T const,1,1> const & x,
            ndarray::Array<T const,1,1> const & y,
            MultiShapeletBasis const & basis,
            MultiShapeletFunction const & psf
        ) {
            _components.reserve(psf.getComponents().size() * basis.getComponentCount());
            for (MultiShapeletBasis::Iterator i = basis.begin(); i != basis.end(); ++i) {
                for (
                    MultiShapeletFunction::ComponentList::const_iterator j = psf.getComponents().begin();
                    j != psf.getComponents().end();
                    ++j
                ) {
                    _components.push_back(
                        std::make_shared< typename RemappedConvolvedShapeletImpl<T>::Factory >(
                            x, y, i->getOrder(), i->getRadius(), i->getMatrix(), *j
                        )
                    );
                }
            }
        }

        virtual int getBasisSize() const { return _components.front()->getBasisSize(); }

        virtual int getDataSize() const { return _components.front()->getDataSize(); }

        virtual int computeWorkspace() const {
            int ws = 0;
            for (FactoryIterator i = _components.begin(); i != _components.end(); ++i) {
                ws = std::max((**i).computeWorkspace(), ws);
            }
            return ws;
        }

        virtual std::shared_ptr<typename MatrixBuilder<T>::Impl> makeBuilderImpl(Workspace & workspace) const {
            Vector builders;
            builders.reserve(_components.size());
            for (FactoryIterator i = _components.begin(); i != _components.end(); ++i) {
                // By copying the workspace here, we prevent the original from having its pointer updated
                // (yet), and hence each component builder starts grabbing workspace arrays from the same
                // point, and they all end up sharing the same space.  That's what we want, because we call
                // them one at a time, and they don't need anything to remain between calls.
                Workspace wsCopy(workspace);
                builders.push_back((**i).makeBuilderImpl(wsCopy));
            }
            // Now, at the end, we increment the workspace by the amount we claimed we needed, which was
            // the maximum needed by any individual components.
            workspace.increment(computeWorkspace());
            return std::make_shared<CompoundImpl>(builders);
        }

    private:
        FactoryVector _components;
    };

    explicit CompoundImpl(Vector const & components) : _components(components) {}

    virtual int getDataSize() const { return _components.front()->getDataSize(); }

    virtual int getBasisSize() const { return _components.front()->getBasisSize(); }

    virtual void buildMatrix(
        ndarray::Array<T,2,-1> const & output,
        afw::geom::ellipses::Ellipse const & ellipse
    ) {
        for (Iterator i = _components.begin(); i != _components.end(); ++i) {
            (**i).buildMatrix(output, ellipse);
        }
    }

private:
    Vector _components;
};

} // anonymous

//===========================================================================================================
//================== MatrixBuilder ==========================================================================
//===========================================================================================================

template <typename T>
MatrixBuilder<T>::MatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order
) {
    *this = MatrixBuilderFactory<T>(x, y, order)();
}

template <typename T>
MatrixBuilder<T>::MatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order,
    ShapeletFunction const & psf
) {
    *this = MatrixBuilderFactory<T>(x, y, order, psf)();
}

template <typename T>
MatrixBuilder<T>::MatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletBasis const & basis
) {
    *this = MatrixBuilderFactory<T>(x, y, basis)();
}

template <typename T>
MatrixBuilder<T>::MatrixBuilder(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletBasis const & basis,
    MultiShapeletFunction const & psf
) {
    *this = MatrixBuilderFactory<T>(x, y, basis, psf)();
}

template <typename T>
int MatrixBuilder<T>::getDataSize() const {
    return _impl->getDataSize();
}

template <typename T>
int MatrixBuilder<T>::getBasisSize() const {
    return _impl->getBasisSize();
}

template <typename T>
ndarray::Array<T,2,-2> MatrixBuilder<T>::allocateOutput() const {
    ndarray::Array<T,2,2> t = ndarray::allocate(getBasisSize(), getDataSize());
    t.deep() = 0.0;
    return t.transpose();
}

template <typename T>
void MatrixBuilder<T>::operator()(
    ndarray::Array<T,2,-1> const & output,
    afw::geom::ellipses::Ellipse const & ellipse
) const {
    _impl->buildMatrix(output, ellipse);
}

template <typename T>
MatrixBuilder<T>::MatrixBuilder(std::shared_ptr<Impl> impl) :
    _impl(impl)
{}


//===========================================================================================================
//================== MatrixBuilderWorkspace =================================================================
//===========================================================================================================

template <typename T>
MatrixBuilderWorkspace<T>::MatrixBuilderWorkspace(int size) {
    std::pair<ndarray::Manager::Ptr,T*> pair = ndarray::SimpleManager<T>::allocate(size);
    _current = pair.second;
    _end = pair.second + size;
    _manager = pair.first;
}


template <typename T>
typename MatrixBuilderWorkspace<T>::Matrix MatrixBuilderWorkspace<T>::makeMatrix(int rows, int cols) {
    Matrix m(_current, rows, cols);
    _current += rows*cols;
    if (_current > _end) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            "Allocated workspace is too small"
        );
    }
    return m;
}

template <typename T>
typename MatrixBuilderWorkspace<T>::Vector MatrixBuilderWorkspace<T>::makeVector(int size) {
    Vector v(_current, size);
    _current += size;
    if (_current > _end) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            "Allocated workspace is too small"
        );
    }
    return v;
}

template <typename T>
void MatrixBuilderWorkspace<T>::increment(int size) {
    _current += size;
}

//===========================================================================================================
//================== MatrixBuilderFactory ===================================================================
//===========================================================================================================

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order
) : _impl(std::make_shared< typename ShapeletImpl<T>::Factory >(x, y, order))
{}

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    int order,
    ShapeletFunction const & psf
) : _impl(std::make_shared< typename ConvolvedShapeletImpl<T>::Factory >(x, y, order, psf))
{}

namespace {

// helper function for the next two ctors: test if a MultiShapeletBasisComponent has unit scaling
// and identity remap matrix
bool isSimple(MultiShapeletBasisComponent const & component) {
    auto matrix = ndarray::asEigenMatrix(component.getMatrix());
    return std::abs(component.getRadius() - 1.0) < std::numeric_limits<double>::epsilon() &&
        matrix.rows() == matrix.cols() &&
        matrix.isApprox(
            Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()),
            std::numeric_limits<double>::epsilon()
        );
}

} // anonymous

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletBasis const & basis
) {
    if (basis.getComponentCount() == 1) {
        MultiShapeletBasisComponent const & component = *basis.begin();
        if (isSimple(component)) {
            _impl = std::make_shared< typename ShapeletImpl<T>::Factory >(x, y, component.getOrder());
        } else {
            _impl = std::make_shared< typename RemappedShapeletImpl<T>::Factory >(
                x, y, component.getOrder(), component.getRadius(), component.getMatrix()
            );
        }
    } else {
        _impl = std::make_shared< typename CompoundImpl<T>::Factory >(x, y, basis);
    }
}

template <typename T>
MatrixBuilderFactory<T>::MatrixBuilderFactory(
    ndarray::Array<T const,1,1> const & x,
    ndarray::Array<T const,1,1> const & y,
    MultiShapeletBasis const & basis,
    MultiShapeletFunction const & psf
) {
    if (basis.getComponentCount() == 1 && psf.getComponents().size() == 1u) {
        ShapeletFunction const & psfComponent = psf.getComponents().front();
        MultiShapeletBasisComponent const & component = *basis.begin();
        if (isSimple(component)) {
            _impl = std::make_shared< typename ConvolvedShapeletImpl<T>::Factory >(
                x, y, component.getOrder(), psfComponent
            );
        } else {
            _impl = std::make_shared< typename RemappedConvolvedShapeletImpl<T>::Factory >(
                x, y, component.getOrder(), component.getRadius(), component.getMatrix(), psfComponent
            );
        }
    } else {
        _impl = std::make_shared< typename CompoundImpl<T>::Factory >(x, y, basis, psf);
    }
}

template <typename T>
int MatrixBuilderFactory<T>::getDataSize() const { return _impl->getDataSize(); }

template <typename T>
int MatrixBuilderFactory<T>::getBasisSize() const { return _impl->getBasisSize(); }

template <typename T>
int MatrixBuilderFactory<T>::computeWorkspace() const { return _impl->computeWorkspace(); }

template <typename T>
MatrixBuilder<T> MatrixBuilderFactory<T>::operator()() const {
    return MatrixBuilder<T>(_impl->makeBuilderImpl());
}

template <typename T>
MatrixBuilder<T> MatrixBuilderFactory<T>::operator()(Workspace & workspace) const {
    return MatrixBuilder<T>(_impl->makeBuilderImpl(workspace));
}

//===========================================================================================================
//================== Explicit Instantiation =================================================================
//===========================================================================================================

#define INSTANTIATE(T)                                          \
    template class MatrixBuilder<T>;                            \
    template class MatrixBuilderFactory<T>;                     \
    template class MatrixBuilderWorkspace<T>

INSTANTIATE(float);
INSTANTIATE(double);

}} // namespace lsst::shapelet
