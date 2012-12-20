#include "lsst/rasmussen/Event.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/geom/Point.h"

namespace lsst {
namespace rasmussen {

Event::Event(afw::image::Image<float> const& im,        // image containing event
             lsst::afw::geom::Point2I const& cen, // central pixel
             int framenum_,                      // frame ID of image
             int chipnum_                        // chip ID for image
            )
{
    if (!im.getBBox(afw::image::PARENT).contains(cen)) {
        throw LSST_EXCEPT(lsst::pex::exceptions::OutOfRangeException,
                          str(boost::format("%d is not contained in %d")
                              % cen % im.getBBox(afw::image::PARENT)));
    }

    framenum = framenum_;
    chipnum = chipnum_;
    x = cen.getX();
    y = cen.getY();
    mode = 0;

    afw::image::Image<float>::xy_locator imData = im.xy_at(cen.getX(), cen.getY());
    int i = 0;
    data[i++] = imData(-1, -1);
    data[i++] = imData( 0, -1);
    data[i++] = imData( 1, -1);
    data[i++] = imData(-1,  0);
    data[i++] = imData( 0,  0);
    data[i++] = imData( 1,  0);
    data[i++] = imData(-1,  1);
    data[i++] = imData( 0,  1);
    data[i++] = imData( 1,  1);
}
        
}
}
