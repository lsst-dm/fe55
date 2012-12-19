#if !defined(LSST_RASMUSSEN_TABLE_H)
#define LSST_RASMUSSEN_TABLE_H
#include "lsst/rasmussen/rv.h"

class HistogramTableBase {
public:
    enum {MAXADU = 4096};
    enum RESET_STYLES { TNONE, T1, T3, T6, };

    HistogramTableBase(const int filter=0x0);
    virtual ~HistogramTableBase() {}
    virtual int process_event(const data_str *ev, int event, int split,
                              RESET_STYLES sty=TNONE, double rst=0.0) = 0;

    void dump_hist(int event, int split, const char *sfile, const char *efile) const;
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
        int *type;
        const int *extr;
        int *hist;
    } table[NMAP];

    int finishEventProcessing(const data_str *ev, const short phe[9], const int map);

    int histo[8][MAXADU];
    int nacc, nnoto;
    unsigned char accmap[NMAP], notomap[NMAP];

    const int _filter;
};

#endif
