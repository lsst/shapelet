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
#ifndef MULTISHAPELET_FitPsf_h_INCLUDED
#define MULTISHAPELET_FitPsf_h_INCLUDED

#include "lsst/ndarray.h"

#include "lsst/afw/math/shapelets.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class FitPsfAlgorithm;

class FitPsfControl : public algorithms::AlgorithmControl {
public:

    LSST_CONTROL_FIELD(innerOrder, int, "Shapelet order of inner expansion (0 == Gaussian)");
    LSST_CONTROL_FIELD(outerOrder, int, "Shapelet order of outer expansion (0 == Gaussian)");

    PTR(FitPsfControl) clone() const { return boost::static_pointer_cast<FitPsfControl>(_clone()); }

    PTR(FitPsfAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)()
    ) const;
    
    FitPsfControl() : algorithms::AlgorithmControl("multishapelet.psf", 2.0), innerOrder(0), outerOrder(0)
        {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief A Multi-Shapelet model for a local PSF.
 *
 *  This is an extension of an elliptical double-Gaussian, where each Gaussian is replaced
 *  with a Gauss-Hermite expansion of arbitrary order (when both orders are zero, it is
 *  exactly a double-Gaussian).  The inner and outer components have the same basis
 *  ellipticity, but the ratio of their radii is not fixed.
 *
 *  See lsst::afw::math::shapelets for the precise definitions of the basis functions.
 */
struct FitPsfModel {

    ndarray::Array<float,1,1> inner; ///< shapelet coefficients of inner expansion
    ndarray::Array<float,1,1> outer; ///< shapelet coefficients of outer expansion
    afw::geom::ellipses::Quadrupole ellipse; ///< ellipse corresponding to inner expansion
    float ratio; ///< radius of outer expansion divided by radius of inner expansion
    bool failed;  ///< set to true if the measurement failed

    /// @brief Construct an uninitialized model with arrays corresponding to the given order.
    FitPsfModel(int innerOrder, int outerOrder);

    /// @brief Construct by extracting saved values from a SourceRecord.
    FitPsfModel(std::string const & name, afw::table::SourceRecord const & source);

    /**
     *  @brief Return a MultiShapeletFunction representation of the model.
     *
     *  The elements will be [inner, outer].
     */
    afw::math::shapelets::MultiShapeletFunction asMultiShapelet(
        afw::geom::Point2D const & center = afw::geom::Point2D()
    ) const;

    template <typename PixelT>
    void evaluate(ndarray::Array<PixelT,2,1> const & array, afw::geom::Point2D const & center) const;

    template <typename PixelT>
    void evaluate(afw::image::Image<PixelT> & image, afw::geom::Point2D const & center) const {
        evaluate(
            image.getArray(),
            afw::geom::Point2D(center.getX() - image.getX0(), center.getY() - image.getY0())
        );
    }

};

class FitPsfAlgorithm : public algorithms::Algorithm {
public:

    /// @brief Construct an algorithm instance and register its fields with a Schema.
    FitPsfAlgorithm(FitPsfControl const & ctrl, afw::table::Schema & schema);

    /// @brief Return the control object
    FitPsfControl const & getControl() const {
        return static_cast<FitPsfControl const &>(algorithms::Algorithm::getControl());
    }

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  This overload accepts the configuration options (inner and outer order) as separate
     *  values, and hence does not require a control object or an Algorithm instance.
     *
     *  @param[in] innerOrder     Shapelet order of the inner expansion (0 == Gaussian)
     *  @param[in] outerOrder     Shapelet order of the outer expansion (0 == Gaussian)
     *  @param[in] image          Postage-stamp image of the PSF.
     *  @param[in] center         Center of the PSF in the image's PARENT coordinate system
     *                            (i.e. xy0 is used).
     */
    template <typename PixelT>
    static FitPsfModel apply(
        int const innerOrder,
        int const outerOrder,
        afw::image::Image<PixelT> const & image,
        afw::geom::Point2D const & center
    );

    /**
     *  @brief Fit a PSF object evaluated at a point.
     *
     *  This overload accepts the configuration options (inner and outer order) as separate
     *  values, and hence does not require a control object or an Algorithm instance.
     *
     *  @param[in] innerOrder     Shapelet order of the inner expansion (0 == Gaussian)
     *  @param[in] outerOrder     Shapelet order of the outer expansion (0 == Gaussian)
     *  @param[in] psf            PSF object
     *  @param[in] center         Point at which to evaluate the PSF.
     */
    static FitPsfModel apply(
        int const innerOrder,
        int const outerOrder,
        afw::detection::Psf const & psf,
        afw::geom::Point2D const & center
    ) {
        return apply(innerOrder, outerOrder, *psf.computeImage(center), center);
    }

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FitPsfAlgorithm);

    afw::table::Key< afw::table::Array<float> > _innerKey;
    afw::table::Key< afw::table::Array<float> > _outerKey;
    afw::table::Key< afw::table::Moments<float> > _ellipseKey;
    afw::table::Key< float > _ratioKey;
    afw::table::Key< afw::table::Flag > _flagKey;
};

inline PTR(FitPsfAlgorithm) FitPsfControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::static_pointer_cast<FitPsfAlgorithm>(_makeAlgorithm(schema, metadata));
}

}}}} // namespace lsst::meas::extensions::multisShapelet

#endif // !MULTISHAPELET_FitPsf_h_INCLUDED
