#if !defined(LSST_RASMUSSEN_EVENT_H)
#define LSST_RASMUSSEN_EVENT_H

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "lsst/rasmussen/rv.h"

namespace lsst {
    namespace afw {
        namespace image {
            template<typename T> class Image;
        }
        namespace geom {
            template<typename T, int N> class Point;
            typedef Point<int, 2> Point2I;
        }
    }
    namespace rasmussen {

        class Event : public data_str {
        public:
            Event(data_str const& ds) : data_str(ds) {}
            Event(lsst::afw::image::Image<float> const& im, ///< image containing event
                  lsst::afw::geom::Point2I const& cen,      ///< central pixel
                  int framenum=-1,                          ///< frame ID of image
                  int chipnum=-1                            ///< chip ID for image
                 );
        };

        std::vector<boost::shared_ptr<Event> > readEventFile(std::string const& fileName);
    }
}
#endif

