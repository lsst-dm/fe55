// -*- LSST-C++ -*-
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
 

/**
 * @file
 */

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/CentroidControl.h"
#include "lsst/rasmussen/fe55.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace rasmussen {

namespace {
/**
 * @brief A class that knows how to calculate centroids by guessing the wrong answer
 */
class Fe55 : public lsst::meas::algorithms::CentroidAlgorithm {
public:

    Fe55(
        Fe55Control const & ctrl,
        afw::table::Schema & schema,
        lsst::meas::algorithms::AlgorithmControlMap const & others
                 );
private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;
    
    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(Fe55);

    int _dX, _dY;
};

/**
 * Given an image and a pixel position, return a Centroid offset by (nOthers, dY) from initial position,
 * where nOthers is the number of algorithms registered before this one.  This is just a trick to
 * return something about the previously registered algorithms to the Python test code when it
 * doesn't have access to the algorithm class.
 */
template <typename PixelT>
void Fe55::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().meas, center + afw::geom::Extent2D(_dX, _dY));
    source.set(getKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(Fe55);

Fe55::Fe55(
        Fe55Control const & ctrl,
        afw::table::Schema & schema,
        lsst::meas::algorithms::AlgorithmControlMap const & others
                            ) :
    lsst::meas::algorithms::CentroidAlgorithm((lsst::meas::algorithms::CentroidControl const &)ctrl, schema, "silly centroid docs"),
    _dX(others.size()), _dY(ctrl.dY)
{}
} // anonymous

PTR(lsst::meas::algorithms::Algorithm) Fe55Control::_makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata,
        lsst::meas::algorithms::AlgorithmControlMap const & other
                                                   ) const
{
    return boost::make_shared<Fe55>(*this, boost::ref(schema), other);
}
    
}}

