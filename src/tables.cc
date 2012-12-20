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
HistogramTableBase::HistogramTableBase(int event, int split,
                                       RESET_STYLES sty, double rst, const int filter) :
    _event(event), _split(split), _filter(filter), _sty(sty), _rst(rst)
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
        t->type = &nsngle;
        t->hist = histo[0];
        t->extr = extra[0];	
    }

    /* load the splus events into table GRADE 1 */
    for (int i = 0; i < sizeof(splus); i++) {
        look_up *t = table + splus[i];
        t->type = &nsplus;
        t->hist = histo[1];
        t->extr = extra[0];	
    }

    /* load the pvert events into table GRADE 2 */
    for (int i = 0; i < sizeof(pvert); i++) {
        look_up *t = table + pvert[i];
        t->type = &npvert;
        t->hist = histo[2];
        t->extr = extra[0];	
    }

    /* load the pleft events into table GRADE 3 */
    for (int i = 0; i < sizeof(pleft); i++) {
        look_up *t = table + pleft[i];
        t->type = &npleft;
        t->hist = histo[3];
        t->extr = extra[0];	
    }

    /* load the prght events into table GRADE 4 */
    for (int i = 0; i < sizeof(prght); i++) {
        look_up *t = table + prght[i];
        t->type = &nprght;
        t->hist = histo[4];
        t->extr = extra[0];	
    }

    /* load the pplus events into table GRADE 5 */
    for (int i = 0; i < sizeof(pplus); i++) {
        look_up *t = table + pplus[i];
        t->type = &npplus;
        t->hist = histo[5];
        t->extr = extra[0];	
    }

    /* load the elnsq events into table GRADE 6 */
    for (int i = 0; i < sizeof(elnsq); i++) {
        look_up *t = table + elnsq[i];
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
    /* make up the array of acceptable maps. */

    for (int i = 0; i < NMAP; i++) accmap[i]=0;
    nacc=0;

    if (filter & 0x01) {
        for (int k=0;k<sizeof(sngle);k++) accmap[nacc+k]=sngle[k];
        nacc+=sizeof(sngle); 
    }
    if (filter & 0x02) {
        for (int k=0;k<sizeof(splus);k++) accmap[nacc+k]=splus[k];
        nacc+=sizeof(splus); 
    }
    if (filter & 0x04) {
        for (int k=0;k<sizeof(pvert);k++) accmap[nacc+k]=pvert[k];
        nacc+=sizeof(pvert); 
    }
    if (filter & 0x08) {
        for (int k=0;k<sizeof(pleft);k++) accmap[nacc+k]=pleft[k];
        nacc+=sizeof(pleft); 
    }
    if (filter & 0x10) {
        for (int k=0;k<sizeof(prght);k++) accmap[nacc+k]=prght[k];
        nacc+=sizeof(prght); 
    }
    if (filter & 0x20) {
        for (int k=0;k<sizeof(pplus);k++) accmap[nacc+k]=pplus[k];
        nacc+=sizeof(pplus); 
    }
    if (filter & 0x40) {
        for (int k=0;k<sizeof(elnsq);k++) accmap[nacc+k]=elnsq[k];
        nacc+=sizeof(elnsq); 
    }
    /* need to make a `not others' array. */
    nnoto=0;
    for (int k=0;k<sizeof(sngle);k++) notomap[nnoto+k]=sngle[k];
    nnoto+=sizeof(sngle);
    for (int k=0;k<sizeof(splus);k++) notomap[nnoto+k]=splus[k];
    nnoto+=sizeof(splus);
    for (int k=0;k<sizeof(pvert);k++) notomap[nnoto+k]=pvert[k];
    nnoto+=sizeof(pvert);
    for (int k=0;k<sizeof(pleft);k++) notomap[nnoto+k]=pleft[k];
    nnoto+=sizeof(pleft);
    for (int k=0;k<sizeof(prght);k++) notomap[nnoto+k]=prght[k];
    nnoto+=sizeof(prght);
    for (int k=0;k<sizeof(pplus);k++) notomap[nnoto+k]=pplus[k];
    nnoto+=sizeof(pplus);
    for (int k=0;k<sizeof(elnsq);k++) notomap[nnoto+k]=elnsq[k];
    nnoto+=sizeof(elnsq);
}

/*
 *  Insert the reset clock correction
 */
void
HistogramTableBase::applyResetClockCorrection(short phe[9])
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

bool
HistogramTableBase::finishEventProcessing(const data_str *ev,
                                          const short phe[9],
                                          const int map)
{
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

    // n.b. grd is in HistogramTableBase
    if (ent->type == &nsngle) {
        grd=0;
    } else if (ent->type == &nsplus) {
        grd=1;
    } else if (ent->type == &npvert) {
        grd=2;
    } else if (ent->type == &npleft) {
        grd=3;
    } else if (ent->type == &nprght) {
        grd=4;
    } else if (ent->type == &npplus) {
        grd=5;
    } else if (ent->type == &nelnsq) {
        grd=6;
    } else if (ent->type == &nother) {
        grd=7;
    } else {
        grd=-1;
    }

    return true;
}

/*
 *  For diagnostic purposes, dump the grade table in CLASSIFY format.
 */
void
HistogramTableBase::dump_table() const
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
HistogramTableBase::dump_head(const char *sfile, int total)
{
    char	line[NAMLEN], *c;
    (void)fprintf(stdout, "!\n");
    (void)fprintf(stdout, "!  QDP Basic Calibration File\n");
    (void)fprintf(stdout, "!\n");
    (void)fprintf(stdout, "!  Working_dir  = %s\n", getcwd((char *)0, NAMLEN));
    (void)fprintf(stdout, "!  Source_file  = %s\n", sfile ? sfile : "unknown");
    (void)fprintf(stdout, "!  Event_thresh = %d\n", _event);
    (void)fprintf(stdout, "!  Split_thresh = %d\n", _split);
    //(void)fprintf(stdout, "!  Total_frames = %d\n", cnt);
    (void)fprintf(stdout, "!  Total_events = %d\n", ntotal);
    (void)fprintf(stdout, "!  Total_pixels = %d\n", (xx - xn)*(yx - yn));
    (void)fprintf(stdout, "!\n");
    (void)fprintf(stdout, "!  Events_below = %d\n", nbevth);
    (void)fprintf(stdout, "!  Events_above = %d\n", noobnd);
    {
        char buff[40];
        if (total < 0) {
            strcpy(buff, "unknown");
        } else {
            sprintf(buff, "%d", total);
        }
                
        (void)fprintf(stdout, "!  Events_input = %s\n", buff);
    }
    (void)fprintf(stdout, "!  PH4_minimum  = %d\n", ev_min);
    (void)fprintf(stdout, "!  PHS_minimum  = %d\n", min_adu);
    (void)fprintf(stdout, "!  PHS_Maximum  = %d\n", max_adu);
    (void)fprintf(stdout, "!  X_Minimum    = %d\n", xn);
    (void)fprintf(stdout, "!  X_Average    = %d\n", xav/(ntotal ? ntotal : 1));
    (void)fprintf(stdout, "!  X_Maximum    = %d\n", xx);
    (void)fprintf(stdout, "!  Y_Minimum    = %d\n", yn);
    (void)fprintf(stdout, "!  Y_Average    = %d\n", yav/(ntotal ? ntotal : 1));
    (void)fprintf(stdout, "!  Y_Maximum    = %d\n", yx);

    (void)fprintf(stdout, "!\n");
    (void)fprintf(stdout, "!  Exclusive grades -- corrected L+Q.\n");
    (void)fprintf(stdout, "!\n");

    if (!sfile || strcmp(sfile, "unknown") == 0) {
        return;
    }
        
    FILE * const fp = fopen(sfile, "r");
    if (fp == NULL) return;

    (void)fprintf(stdout, "!\n");
    (void)fprintf(stdout, "!  Experimental parameters\n");
    (void)fprintf(stdout, "!\n");
    while (fgets(line, NAMLEN-1, fp)) {
        if (line[0] == '#') continue;
        (void)fprintf(stdout, "!  %s", line);
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
HistogramTableBase::dump_hist(const char *sfile) const
{
    (void)printf("!\n");
    (void)printf("!  QDP Header follows\n");
    (void)printf("!\n");
    (void)printf("lab top Event = %d Split = %d Source = %s\n",
                 _event, _split, sfile ? sfile : "unknown");
    if (_efile[0]) (void)printf("lab file %s", _efile);
    (void)printf("lab g1 Pulse Height (ADU)\n");
    (void)printf("lab rot\n");
    (void)printf("lab g2 N(S)\nlab g3 N(S+)\n");
    (void)printf("lab g4 N(Pv)\nlab g5 N(Pl)\n");
    (void)printf("lab g6 N(Pr)\nlab g7 N(P+)\n");
    (void)printf("lab g8 N(L+Q)\nlab g9 N(O)\n");
    (void)printf("csize 0.75\n");

    const int EXTADU = 8;
    const int tmp_min_2ct = (min_2ct <        EXTADU) ?      0 : min_2ct - EXTADU;
    const int tmp_max_2ct = (max_2ct >= MAXADU - EXTADU) ? MAXADU : max_2ct + 1 + EXTADU;
    (void)printf("res x %d %d\n", tmp_min_2ct, tmp_max_2ct);
    (void)printf("res y2 1\nres y3 1\n");
    (void)printf("res y4 1\nres y5 1\n");
    (void)printf("res y6 1\nres y7 1\n");
    (void)printf("res y8 1\nres y9 1\n");

    (void)printf("error y sq 2 3 4 5 6 7 8 9\n");
    (void)printf("log y on\n");
    (void)printf("plot vert\n");
    (void)printf("!\n");
    (void)printf("!  Histogram data follows\n");
    (void)printf("!\n");
    (void)printf("!  PHA\tS\tS+\tPv\tPl\tPr\tP+\tL+Q\tO\n");
    (void)printf("!  TOT\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                  nsngle,nsplus,npvert,npleft,nprght,npplus,nelnsq,nother);
    (void)printf("!\n");

    int tmp_min_adu = (min_adu <         EXTADU) ?      0 : min_adu - EXTADU;
    int tmp_max_adu = (max_adu >= MAXADU-EXTADU) ? MAXADU : max_adu + 1 + EXTADU;
    for (int i = tmp_min_adu; i < tmp_max_adu; i++) {
        (void)printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i,
                     histo[0][i], histo[1][i], histo[2][i], histo[3][i],
                     histo[4][i], histo[5][i], histo[6][i], histo[7][i]);
    }
}
