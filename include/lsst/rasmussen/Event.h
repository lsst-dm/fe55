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
            // Event's "Grade"
            enum Grade { UNKNOWN=-1,    ///< Unknown
                         SINGLE=0,      ///< Single pixel
                         SINGLE_P_CORNER=1, ///< Single pixel + corner(s)
                         VERTICAL_SPLIT=2,  ///< Vertical split (+ detached corner(s))
                         LEFT_SPLIT=3,      ///< Left split (+ detached corner(s))
                         RIGHT_SPLIT=4,     ///< Right split (+ detached corner(s))
                         SINGLE_SIDED_P_CORNER=5, ///< Single-sided split +  touched corner
                         ELL_SQUARE_P_CORNER=6,   ///< L or square (+ detached corner)
                         OTHER=7                  ///< all others
            };

            Event(data_str const& ds) : data_str(ds), grade(UNKNOWN), sum(0.0), p9(0.0) {}
            Event(lsst::afw::image::Image<float> const& im, ///< image containing event
                  lsst::afw::geom::Point2I const& cen,      ///< central pixel
                  int framenum=-1,                          ///< frame ID of image
                  int chipnum=-1                            ///< chip ID for image
                 );

            Grade grade;                  ///< Event's grade
            float sum;                    ///< Sum of counts in Event
            float p9;                     ///< Event's "P9" sum
        };

        std::vector<boost::shared_ptr<Event> > readEventFile(std::string const& fileName);
    }
}
#endif

