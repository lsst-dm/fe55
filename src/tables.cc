/*
 *  rv_ev2pcf.c -- Convert evlist to basic calibration file
 *
 *	Author:		gbc@space.mit.edu	17 Oct 1992
 *
 *	The primary calibration file is a QDP file giving a 
 *	histogram of the number of events of each Astro-D
 *	grade, preceeded by QDP/PLT plotting commands.  A
 *	header giving the experimental source parameters
 *	is also included in the form of QDP comment lines.
 *
 *	Modified for new SIS grades	gbc	02 Dec 1992
 *	Modified for Reset correction	gbc	19 Mar 1993
 */
#include <limits>
#include <algorithm>
#include "lsst/rasmussen/tables.h"

/*
 *  Initialize the histogram tables.  Each entry of the table needs
 *  to know which counter to increment, which histogram to increment
 *  and which extra pixels should be included in the summed pha,
 *  which only occurs for the L, Q and Other grades.
 */
HistogramTable::HistogramTable(int event, int split,
                                       RESET_STYLES sty, double rst, const int filter,
                                       calctype do_what) :
    histo(ndarray::allocate(ndarray::makeVector(8, MAXADU))),
    _event(event), _split(split), _filter(filter), _do_what(do_what), _efile(""), _sty(sty), _rst(rst)
{
    static
    const int extra[][4] = {  {4,4,4,4},
                              {0,4,4,4}, {2,4,4,4}, {6,4,4,4}, {8,4,4,4},
                              {0,2,4,4}, {0,6,4,4}, {6,8,4,4}, {8,2,4,4},
                              {0,2,6,8} };
    static
    const unsigned char emask[] = { 0x00,
                                    0x0a,0x12,0x48,0x50,
                                    0x1a,0x4a,0x58,0x52,
                                    0x5a };
    const unsigned char sngle[] = { 0x00 };					/* 0 */
    const unsigned char splus[] = { 0x01,0x04,0x20,0x80,0x05,0x21,0x81,		/* 1 */
				    0x24,0x84,0xa0,0x25,0x85,0xa4,0xa1,0xa5 };
    const unsigned char pvert[] = { 0x02,0x40,0x22,0x82,0x41,0x44,0x45,0xa2 };	/* 2 */
    const unsigned char pleft[] = { 0x08,0x0c,0x88,0x8c };			/* 3 */
    const unsigned char prght[] = { 0x10,0x30,0x11,0x31 };			/* 4 */
/*
 *  const unsigned char phorz[] = { 0x08,0x10,0x0c,0x88,0x30,0x11,0x8c,0x31 };
 */
    const unsigned char pplus[] = { 0x03,0x06,0x09,0x28,0x60,0xc0,0x90,0x14,	/* 5 */
				    0x83,0x26,0x89,0x2c,0x64,0xc1,0x91,0x34,
				    0x23,0x86,0x0d,0xa8,0x61,0xc4,0xb0,0x15,
				    0xa3,0xa6,0x8d,0xac,0x65,0xc5,0xb1,0x35 };
/*
 *  const unsigned char elish[] = { 0x12,0x32,0x50,0x51,0x48,0x4c,0x0a,0x8a };
 *  const unsigned char squar[] = { 0x16,0xd0,0x68,0x0b,0x36,0xd1,0x6c,0x8b };
 */
    const unsigned char elnsq[] = { 0x12,0x32,0x50,0x51,0x48,0x4c,0x0a,0x8a,	/* 6 */
				    0x16,0xd0,0x68,0x0b,0x36,0xd1,0x6c,0x8b };

#if 0							  // not actually needed as an array
    const unsigned char other[NMAP] = { all of the rest };			/* 7 */
#endif

    /* initialize everything in sight */
    nsngle = nsplus = npvert = npleft = nprght = npplus =
        nelnsq = nother = ntotal = noobnd = nbevth = 0;

    ev_min = MAXADU; xav = 0; yav = 0;
    min_adu = MAXADU; max_adu = 0;
    min_2ct = MAXADU; max_2ct = 0;
    xn = yn = std::numeric_limits<int>::max();
    xx = yx = 0;

    std::fill(&histo[0][0], &histo[0][0] + 8*MAXADU, 0);

    /* load the sngle events into table GRADE 0 */
    for (int i = 0; i < sizeof(sngle); i++) {
        look_up *t = table + sngle[i];
        t->grade = lsst::rasmussen::Event::SINGLE;
        t->type = &nsngle;
        t->hist = histo[0];
        t->extr = extra[0];	
    }

    /* load the splus events into table GRADE 1 */
    for (int i = 0; i < sizeof(splus); i++) {
        look_up *t = table + splus[i];
        t->grade = lsst::rasmussen::Event::SINGLE_P_CORNER;
        t->type = &nsplus;
        t->hist = histo[1];
        t->extr = extra[0];	
    }

    /* load the pvert events into table GRADE 2 */
    for (int i = 0; i < sizeof(pvert); i++) {
        look_up *t = table + pvert[i];
        t->grade = lsst::rasmussen::Event::VERTICAL_SPLIT;
        t->type = &npvert;
        t->hist = histo[2];
        t->extr = extra[0];	
    }

    /* load the pleft events into table GRADE 3 */
    for (int i = 0; i < sizeof(pleft); i++) {
        look_up *t = table + pleft[i];
        t->grade = lsst::rasmussen::Event::LEFT_SPLIT;
        t->type = &npleft;
        t->hist = histo[3];
        t->extr = extra[0];	
    }

    /* load the prght events into table GRADE 4 */
    for (int i = 0; i < sizeof(prght); i++) {
        look_up *t = table + prght[i];
        t->grade = lsst::rasmussen::Event::RIGHT_SPLIT;
        t->type = &nprght;
        t->hist = histo[4];
        t->extr = extra[0];	
    }

    /* load the pplus events into table GRADE 5 */
    for (int i = 0; i < sizeof(pplus); i++) {
        look_up *t = table + pplus[i];
        t->grade = lsst::rasmussen::Event::SINGLE_SIDED_P_CORNER;
        t->type = &npplus;
        t->hist = histo[5];
        t->extr = extra[0];	
    }

    /* load the elnsq events into table GRADE 6 */
    for (int i = 0; i < sizeof(elnsq); i++) {
        look_up *t = table + elnsq[i];
        t->grade = lsst::rasmussen::Event::ELL_SQUARE_P_CORNER;
        t->type = &nelnsq;
        t->hist = histo[6];
        t->extr = extra[0];	
        for (int b = (0x5a & elnsq[i]), j = 0; j < sizeof(emask); j++)
            if (b == emask[j]) {
                t->extr = extra[j];
                break;
            }
    }

    /* load the other events into table GRADE 7 */
    for (int i = 0; i < NMAP; i++) {
        look_up *t = table + i;
        if (t->type) continue;		/* already loaded */
        t->grade = lsst::rasmussen::Event::OTHER;
        t->type = &nother;
        t->hist = histo[7];
        t->extr = extra[0];	
        /*
         *  In this version, included corners are
         *  ignored in all of the grade 7 events.
         *
         for (b = 0x5a & (unsigned char)i, j = 0;
         j < sizeof(emask); j++)
         if (b == emask[j]) {
         t->extr = extra[j];
         break;
         }
        */
    }
}

const int HistogramTable::MAXADU = 4096;

/*
 *  Insert the reset clock correction
 */
void
HistogramTable::applyResetClockCorrection(short phe[9])
{
    switch (_sty) {
      case T6:
        phe[7] -= phe[6]*_rst;
        phe[4] -= phe[3]*_rst;
        phe[1] -= phe[0]*_rst; // fall through
      case T3:
        phe[8] -= phe[7]*_rst;
        phe[2] -= phe[1]*_rst; // fall through
      case T1:
        phe[5] -= phe[4]*_rst;
        break;
      case TNONE:
        break;
    }
}

/*********************************************************************************************************/
/*
 *  Accumulate the num events in the tables
 */
bool
HistogramTable::process_event(lsst::rasmussen::Event *ev
                                 )
{
    /*
     *  Get some gross event parameters
     */
    if (ev->data[4] < ev_min) ev_min = ev->data[4];
    if (ev->data[4] < _event) {
        nbevth++;
        ev->grade = lsst::rasmussen::Event::UNKNOWN; // We don't know map yet.
        return false;
    }

    short phe[9];
    std::copy(ev->data, ev->data + 9, phe);

    applyResetClockCorrection(phe);
    /*
     *  Characterize event & accumulate most of pha
     */
    unsigned char map = 0;

    ev->p9 = 0;
    sum = 0;                            // in HistogramTable
    for (int j = 0; j < 9; j++) {
        const short phj = phe[j];

        switch (_do_what) {
          case P_9:
            ev->p9 += phj;
            break;
          case P_1357:
            if (j == 1 || j == 3 || j == 5 || j == 7 || j == 4) ev->p9 += phj;
            break;
          case P_17:
            if (j == 1 || j == 7 || j == 4) ev->p9 += phj;
            break;
          case P_35:
            if (j == 3 || j == 5 || j == 4) ev->p9 += phj;
            break;
          case P_LIST:
            break;
        }

        if (phj < _split && j != 4) {
            phe[j] = 0;
            continue;
        }
        switch (j) {
          case 0: map |= 0x01;           ; break;
          case 1: map |= 0x02; sum += phj; break;
          case 2: map |= 0x04;           ; break;
          case 3: map |= 0x08; sum += phj; break;
          case 4: 	       sum += phj; break;
          case 5: map |= 0x10; sum += phj; break;
          case 6: map |= 0x20;           ; break;
          case 7: map |= 0x40; sum += phj; break;
          case 8: map |= 0x80;           ; break;
        }
    }
    ev->grade = table[map].grade;
    /* 
     *  grade is identified. check with _filter to see whether  to pass it on or not.
     */
    if (ev->grade == lsst::rasmussen::Event::UNKNOWN ||
        ((1 << static_cast<int>(ev->grade)) & _filter) == 0x0) {
        return false;
    }

    /*
     *  Finish pha with extra pixels of L, Q, and O events
     */
    look_up *const ent = &table[map];
    const int *xtr = ent->extr;
    for (int j = 0; xtr[j] != 4 && j < 4; j++) sum += phe[xtr[j]];
    /*
     *  Accumulate statistics and various bounds
     */
    if (sum >= MAXADU) { noobnd++;  return false; }
    if (sum > max_adu) max_adu = sum;
    if (sum < min_adu) min_adu = sum;
    if (ev->x < xn) xn = ev->x;
    if (ev->x > xx) xx = ev->x;
    if (ev->y < yn) yn = ev->y;
    if (ev->y > yx) yx = ev->y;
    xav += ev->x;
    yav += ev->y;
    ntotal += 1;
    *ent->type += 1;
    const int hsum = ent->hist[sum] += 1;
    if (hsum > 2) {
        if (sum > max_2ct) max_2ct = sum;
        if (sum < min_2ct) min_2ct = sum;
    }

    return true;
}

/*
 *  For diagnostic purposes, dump the grade table in CLASSIFY format.
 */
void
HistogramTable::dump_table() const
{
    for (int i = 0, j; i != NMAP; ++i) {
        if      ( table[i].type == &nsngle ) j = 0;
	else if ( table[i].type == &nsplus ) j = 1;
        else if ( table[i].type == &npvert ) j = 2;
        else if ( table[i].type == &npleft ) j = 3;
        else if ( table[i].type == &nprght ) j = 4;
        else if ( table[i].type == &npplus ) j = 5;
        else if ( table[i].type == &nelnsq ) j = 6;
        else if ( table[i].type == &nother ) j = 7;
        else                                 j = 9;

        (void)fprintf(stderr, "%d,", j);
        if (i%16 == 15) (void)fprintf(stderr, "\n");
    }
}


/*
 *  Dump the basic calibration file header
 */
void
HistogramTable::dump_head(FILE *fd,
                              const char *sfile, int total)
{
    char	line[NAMLEN], *c;
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  QDP Basic Calibration File\n");
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  Working_dir  = %s\n", getcwd((char *)0, NAMLEN));
    (void)fprintf(fd, "!  Source_file  = %s\n", sfile ? sfile : "unknown");
    (void)fprintf(fd, "!  Event_thresh = %d\n", _event);
    (void)fprintf(fd, "!  Split_thresh = %d\n", _split);
    //(void)fprintf(fd, "!  Total_frames = %d\n", cnt);
    (void)fprintf(fd, "!  Total_events = %d\n", ntotal);
    (void)fprintf(fd, "!  Total_pixels = %d\n", (xx - xn)*(yx - yn));
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  Events_below = %d\n", nbevth);
    (void)fprintf(fd, "!  Events_above = %d\n", noobnd);
    {
        char buff[40];
        if (total < 0) {
            strcpy(buff, "unknown");
        } else {
            sprintf(buff, "%d", total);
        }
                
        (void)fprintf(fd, "!  Events_input = %s\n", buff);
    }
    (void)fprintf(fd, "!  PH4_minimum  = %d\n", ev_min);
    (void)fprintf(fd, "!  PHS_minimum  = %d\n", min_adu);
    (void)fprintf(fd, "!  PHS_Maximum  = %d\n", max_adu);
    (void)fprintf(fd, "!  X_Minimum    = %d\n", xn);
    (void)fprintf(fd, "!  X_Average    = %d\n", xav/(ntotal ? ntotal : 1));
    (void)fprintf(fd, "!  X_Maximum    = %d\n", xx);
    (void)fprintf(fd, "!  Y_Minimum    = %d\n", yn);
    (void)fprintf(fd, "!  Y_Average    = %d\n", yav/(ntotal ? ntotal : 1));
    (void)fprintf(fd, "!  Y_Maximum    = %d\n", yx);

    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  Exclusive grades -- corrected L+Q.\n");
    (void)fprintf(fd, "!\n");

    if (!sfile || strcmp(sfile, "unknown") == 0) {
        return;
    }
        
    FILE * const fp = fopen(sfile, "r");
    if (fp == NULL) return;

    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  Experimental parameters\n");
    (void)fprintf(fd, "!\n");
    while (fgets(line, NAMLEN-1, fp)) {
        if (line[0] == '#') continue;
        (void)fprintf(fd, "!  %s", line);
        if (!strncmp(line, "evlist  = ", 10)) {
            (void)strcpy(_efile, &line[10]);
            for (c = _efile; *c; c++)
                if (*c == '\t') *c = ' '; 
        }
    }
    (void)fclose(fp);
}

/*
 *  Dump the QDP header histogram table
 */
void
HistogramTable::dump_hist(FILE *fd,
                              const char *sfile) const
{
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  QDP Header follows\n");
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "lab top Event = %d Split = %d Source = %s\n",
                 _event, _split, sfile ? sfile : "unknown");
    if (_efile[0]) (void)fprintf(fd, "lab file %s", _efile);
    (void)fprintf(fd, "lab g1 Pulse Height (ADU)\n");
    (void)fprintf(fd, "lab rot\n");
    (void)fprintf(fd, "lab g2 N(S)\nlab g3 N(S+)\n");
    (void)fprintf(fd, "lab g4 N(Pv)\nlab g5 N(Pl)\n");
    (void)fprintf(fd, "lab g6 N(Pr)\nlab g7 N(P+)\n");
    (void)fprintf(fd, "lab g8 N(L+Q)\nlab g9 N(O)\n");
    (void)fprintf(fd, "csize 0.75\n");

    const int EXTADU = 8;
    const int tmp_min_2ct = (min_2ct <        EXTADU) ?      0 : min_2ct - EXTADU;
    const int tmp_max_2ct = (max_2ct >= MAXADU - EXTADU) ? MAXADU : max_2ct + 1 + EXTADU;
    (void)fprintf(fd, "res x %d %d\n", tmp_min_2ct, tmp_max_2ct);
    (void)fprintf(fd, "res y2 1\nres y3 1\n");
    (void)fprintf(fd, "res y4 1\nres y5 1\n");
    (void)fprintf(fd, "res y6 1\nres y7 1\n");
    (void)fprintf(fd, "res y8 1\nres y9 1\n");

    (void)fprintf(fd, "error y sq 2 3 4 5 6 7 8 9\n");
    (void)fprintf(fd, "log y on\n");
    (void)fprintf(fd, "plot vert\n");
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  Histogram data follows\n");
    (void)fprintf(fd, "!\n");
    (void)fprintf(fd, "!  PHA\tS\tS+\tPv\tPl\tPr\tP+\tL+Q\tO\n");
    (void)fprintf(fd, "!  TOT\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                  nsngle,nsplus,npvert,npleft,nprght,npplus,nelnsq,nother);
    (void)fprintf(fd, "!\n");

    int tmp_min_adu = (min_adu <         EXTADU) ?      0 : min_adu - EXTADU;
    int tmp_max_adu = (max_adu >= MAXADU-EXTADU) ? MAXADU : max_adu + 1 + EXTADU;
    for (int i = tmp_min_adu; i < tmp_max_adu; i++) {
        (void)fprintf(fd, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i,
                      histo[0][i], histo[1][i], histo[2][i], histo[3][i],
                      histo[4][i], histo[5][i], histo[6][i], histo[7][i]);
    }
}
