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

#include "lsst/meas/extensions/multiShapelet/FitPsf.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

namespace {

int computeOrder(int size2d) {
    // Could just use f.p. square roots here, but we only care about small numbers
    // and I don't think about how to be robust against round-off error issues.
    for (int o=0, s=1; s < size2d; s += (++o + 1)) {
        if (s == size2d) return o;
    }
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        (boost::format("%d is not a valid shapelet vector size") % size2d).str()
    );
}

template <typename PixelT>
void fillMultiShapeletImage(
    ndarray::Array<PixelT,2,1> const & array, 
    afw::math::shapelets::MultiShapeletFunction const & msf
) {
    afw::math::shapelets::MultiShapeletFunctionEvaluator ev(msf);
    afw::geom::Point2D point;
    for (
        typename ndarray::Array<PixelT,2,1>::Iterator rowIter = array.begin();
        rowIter != array.end();
        ++rowIter, point.setY(point.getY() + 1.0) // FIXME: overload Point::getY to return reference
    ) {
        for (
            typename ndarray::Array<PixelT,2,1>::Reference::Iterator pixIter = rowIter->begin();
            pixIter != rowIter->end();
            ++pixIter, point.setX(point.getX() + 1.0) // FIXME: overload Point::getX to return reference
        ) {
            *pixIter = ev(point);
        }
        point.setY(0.0);
    }
}

} // anonymous

FitPsfModel::FitPsfModel(FitPsfControl const & ctrl) :
    inner(ndarray::allocate(afw::math::shapelets::computeSize(ctrl.innerOrder))),
    outer(ndarray::allocate(afw::math::shapelets::computeSize(ctrl.outerOrder))),
    ellipse(),
    radiusRatio(ctrl.radiusRatio),
    failed(false)
{
    inner.deep() = 0.0;
    outer.deep() = 0.0;
}

FitPsfModel::FitPsfModel(FitPsfControl const & ctrl, afw::table::SourceRecord const & source) :
    inner(ndarray::allocate(afw::math::shapelets::computeSize(ctrl.innerOrder))),
    outer(ndarray::allocate(afw::math::shapelets::computeSize(ctrl.outerOrder))),
    ellipse(),
    radiusRatio(ctrl.radiusRatio),
    failed(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    inner.deep() = source.get(s.find< afw::table::Array<float> >("inner").key);
    outer.deep() = source.get(s.find< afw::table::Array<float> >("outer").key);
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failed = source.get(s.find<afw::table::Flag>("flags").key);
}
  

afw::math::shapelets::MultiShapeletFunction FitPsfModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    afw::math::shapelets::MultiShapeletFunction::ElementList elements;
    afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
    // FIXME: should have data-type cast functionality in ndarray
    ndarray::Array<afw::math::shapelets::Pixel,1,1> sInner(ndarray::allocate(inner.getSize<0>()));
    ndarray::Array<afw::math::shapelets::Pixel,1,1> sOuter(ndarray::allocate(outer.getSize<0>()));
    sInner.deep() = inner;
    sOuter.deep() = outer;
    elements.push_back(
        afw::math::shapelets::ShapeletFunction(
            computeOrder(inner.getSize<0>()),
            afw::math::shapelets::HERMITE,
            fullEllipse,
            sInner
        )
    );
    fullEllipse.scale(radiusRatio);
    elements.push_back(
        afw::math::shapelets::ShapeletFunction(
            computeOrder(outer.getSize<0>()),
            afw::math::shapelets::HERMITE,
            fullEllipse,
            sOuter
        )
    );
    return afw::math::shapelets::MultiShapeletFunction(elements);
}

template <typename PixelT>
void FitPsfModel::evaluate(
    ndarray::Array<PixelT,2,1> const & array, afw::geom::Point2D const & center
) const {
    afw::math::shapelets::MultiShapeletFunction msf = asMultiShapelet(center);
    fillMultiShapeletImage(array, msf);
}

FitPsfAlgorithm::FitPsfAlgorithm(FitPsfControl const & ctrl, afw::table::Schema & schema) :
    algorithms::Algorithm(ctrl),
    _innerKey(schema.addField< afw::table::Array<float> >(
                  ctrl.name + ".inner",
                  "Gauss-Hermite coefficients of the inner expansion (see afw::math::shapelets)",
                  afw::math::shapelets::computeSize(ctrl.innerOrder)
              )),
    _outerKey(schema.addField< afw::table::Array<float> >(
                  ctrl.name + ".outer",
                  "Gauss-Hermite coefficients of the outer expansion (see afw::math::shapelets)",
                  afw::math::shapelets::computeSize(ctrl.outerOrder)
              )),
    _ellipseKey(schema.addField< afw::table::Moments<float> >(
                    ctrl.name + ".ellipse",
                    "Ellipse corresponding to the inner expansion"
                )),
    _flagKey(schema.addField<afw::table::Flag>(
                 ctrl.name + ".flags",
                 "set if the multi-shapelet PSF fit was unsuccessful"
             ))
{}

template <typename PixelT>
void FitPsfAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitPsfAlgorithm without a PSF."
        );
    }
    source.set(_flagKey, true);
    FitPsfModel model = apply(getControl(), *exposure.getPsf(), center);
    source[_innerKey] = model.inner;
    source[_outerKey] = model.outer;
    source.set(_ellipseKey, model.ellipse);
    source.set(_flagKey, model.failed);
}

template <typename PixelT>
FitPsfModel FitPsfAlgorithm::apply(
    FitPsfControl const & ctrl,
    afw::image::Image<PixelT> const & image,
    afw::geom::Point2D const & center
) {
    // TODO
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitPsfAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
