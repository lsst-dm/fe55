#if !defined(LSST_RASMUSSEN_TABLE_H)
#define LSST_RASMUSSEN_TABLE_H
#include "lsst/rasmussen/rv.h"

class HistogramTableBase {
public:
    enum {MAXADU = 4096};
    enum RESET_STYLES { TNONE, T1, T3, T6, };

    HistogramTableBase(const int filter=0x0);
    void dump_table() const;

    int		nsngle,nsplus,npvert,npleft,nprght,npplus,
		nelnsq,nother,ntotal,noobnd,nbevth;
    int		ev_min, xav, yav;
    short	min_adu, max_adu;
    short	min_2ct, max_2ct;
    short	xn, xx, yn, yx;
protected:
    enum { NMAP = 256 };
    struct look_up {
        int *type;
        const int *extr;
        int *hist;
    } table[NMAP];

    int histo[8][MAXADU];
    int nacc, nnoto;
    unsigned char accmap[NMAP], notomap[NMAP];

    const int _filter;
};

#endif
