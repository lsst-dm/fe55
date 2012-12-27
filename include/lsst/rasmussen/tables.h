#if !defined(LSST_RASMUSSEN_TABLE_H)
#define LSST_RASMUSSEN_TABLE_H
#include "ndarray.h"
#include "lsst/rasmussen/Event.h"

/**
 * \brief Hello World
 */
class HistogramTableBase {
public:
    static const int MAXADU;
    enum RESET_STYLES { TNONE, T1, T3, T6, };

    HistogramTableBase(int event=0, int split=0,
                       RESET_STYLES sty=TNONE, double rst=0.0, const int filter=0x0);
    virtual ~HistogramTableBase() {}
    virtual bool process_event(lsst::rasmussen::Event *ev) { return false; }

    void dump_head(FILE *fd=stdout, const char *sfile=NULL, int total=-1);
    void dump_hist(FILE *fd=stdout, const char *sfile=NULL) const;
    void dump_table() const;

    int		nsngle,nsplus,npvert,npleft,nprght,npplus,
		nelnsq,nother,ntotal,noobnd,nbevth;
    int		ev_min, xav, yav;
    int		min_adu, max_adu;
    int		min_2ct, max_2ct;
    int		xn, xx, yn, yx;

#define USE_NDARRAY 1
#if USE_NDARRAY
    ndarray::Array<int, 2, 2> histo;
#else
    int histo[8][MAXADU];
#endif

    // Values set by process_event
    int sum;                            // should be float?  But it's used as an array index
protected:
    enum { NMAP = 256 };
    struct look_up {
        look_up() : type(0), extr(0),
#if USE_NDARRAY
                    hist()
#else
                    hist(0)
#endif
        {}

        int *type;
        const int *extr;
#if USE_NDARRAY
    ndarray::Array<int, 1, 1> hist;
#else
        int *hist;
#endif
    } table[NMAP];

    void applyResetClockCorrection(short phe[9]);
    bool finishEventProcessing(lsst::rasmussen::Event *ev, const short phe[9], const int map);
    lsst::rasmussen::Event::Grade setGrdFromType(const int map);

    int nacc, nnoto;
    unsigned char accmap[NMAP], notomap[NMAP];

    int _event;
    int _split;
    const int _filter;
private:
    enum { NAMLEN = 512 };

    char _efile[NAMLEN];                // name of the electronics param file, found in the sfile.  Ughh
    RESET_STYLES _sty;
    double _rst;
};

/********************************************************************************************************/

class HistogramTableGflt : public HistogramTableBase {
public:
    HistogramTableGflt(const int filter=0x0, int event=0, int split=0,
                       RESET_STYLES sty=TNONE, double rst=0.0) :
        HistogramTableBase(event, split, sty, rst, filter) {}
    virtual bool process_event(lsst::rasmussen::Event *ev);
};

/*********************************************************************************************************/

class HistogramTableXygpx : public HistogramTableBase {
public:
    enum calctype { P_9,
                    P_17,
                    P_35,
                    P_1357,
                    P_LIST,             // for the "total"
    };

    HistogramTableXygpx(calctype do_what=P_LIST, int event=0, int split=0,
                        RESET_STYLES sty=TNONE, double rst=0.0) :
        HistogramTableBase(event, split, sty, rst), _do_what(do_what) {}

    virtual bool process_event(lsst::rasmussen::Event *ev);
private:
    const calctype _do_what;
};

#endif
