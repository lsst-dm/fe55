#include <cstdio>
#include "lsst/rasmussen/Event.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/geom/Point.h"

namespace lsst {
namespace rasmussen {

Event::Event(afw::image::Image<float> const& im,        // image containing event
             lsst::afw::geom::Point2I const& cen,       // central pixel
             int framenum_,                             // frame ID of image
             int chipnum_                               // chip ID for image
            ) : grade(UNKNOWN), sum(0.0), p9(0.0)

{
    if (!im.getBBox(afw::image::PARENT).contains(cen - afw::geom::ExtentI(1, 1)) ||
        !im.getBBox(afw::image::PARENT).contains(cen + afw::geom::ExtentI(1, 1))) {
        throw LSST_EXCEPT(lsst::pex::exceptions::OutOfRangeException,
                          str(boost::format("%d is too close to the edge of image of size %d")
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

std::vector<PTR(Event)>
readEventFile(std::string const& fileName)
{
    FILE *fp = fopen(fileName.c_str(), "r");
    if (!fp) {
        throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException,
                          str(boost::format("Unable to open %s for read")
                              % fileName));
    }

    std::vector<PTR(Event)> events;

    data_str event;
    while (fread((void *)&event, sizeof(data_str), 1, fp) > 0) {
        events.push_back(PTR(Event)(new Event(event)));
    }

    fclose(fp);

    return events;
}

}}
