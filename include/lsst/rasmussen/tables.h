#if !defined(LSST_RASMUSSEN_TABLE_H)
#define LSST_RASMUSSEN_TABLE_H
#include "lsst/rasmussen/rv.h"

class HistogramTableBase {
public:
    enum {MAXADU = 4096};
    enum RESET_STYLES { TNONE, T1, T3, T6, };

    HistogramTableBase(int event=0, int split=0,
                       RESET_STYLES sty=TNONE, double rst=0.0, const int filter=0x0);
    virtual ~HistogramTableBase() {}
    virtual bool process_event(const data_str *ev) { return false; }

    void dump_head(const char *sfile=NULL, int total=-1);
    void dump_hist(const char *sfile=NULL) const;
    void dump_table() const;

    int		nsngle,nsplus,npvert,npleft,nprght,npplus,
		nelnsq,nother,ntotal,noobnd,nbevth;
    int		ev_min, xav, yav;
    int		min_adu, max_adu;
    int		min_2ct, max_2ct;
    int		xn, xx, yn, yx;
    // Values set by process_event
    int grd;                            // the event's grade
    int sum;                            // should be float?  But it's used as an array index
protected:
    enum { NMAP = 256 };
    struct look_up {
        look_up() : type(0), extr(0), hist(0) {}

        int *type;
        const int *extr;
        int *hist;
    } table[NMAP];

    void applyResetClockCorrection(short phe[9]);
    bool finishEventProcessing(const data_str *ev, const short phe[9], const int map);

    int histo[8][MAXADU];
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
    virtual bool process_event(const data_str *ev);
};

/*********************************************************************************************************/

class HistogramTableXygpx : public HistogramTableBase {
public:
    enum calctype { p_9,
                    p_17,
                    p_35,
                    p_1357,
                    p_list,             // for the "total"
    };

    HistogramTableXygpx(calctype do_what=p_list, int event=0, int split=0,
                        RESET_STYLES sty=TNONE, double rst=0.0) :
        HistogramTableBase(event, split, sty, rst), _do_what(do_what) {}

    virtual bool process_event(const data_str *ev);

    // Value set by process_event
    int p9;

private:
    const calctype _do_what;
};

#endif
