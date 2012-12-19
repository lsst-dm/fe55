/*
 *  rv_gradefilt.c -- filter evlist based on ASUKA grades
 *
 *    orig. Author:   gbc@space.mit.edu       17 Oct 1992
 *
 *	The primary calibration file is a QDP file giving a 
 *	histogram of the number of events of each Astro-D
 *	grade, preceeded by QDP/PLT plotting commands.  A
 *	header giving the experimental source parameters
 *	is also included in the form of QDP comment lines.
 *
 *	Modified for new SIS grades	gbc	02 Dec 1992
 *	Modified for Reset correction	gbc	19 Mar 1993
 *
 *      Hacked into filter from rv_ev2pcf by Andy Rasmussen, 21 MAR 1993
 *           made to pass named grades on the standard output taken from
 *           the standard inut.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "lsst/rasmussen/rv.h"

namespace {
    const int MAXADU = 4096;

    struct look_up {
        int *type;
        const int *extr;
        int *hist;
    } table[256];

    int		nsngle,nsplus,npvert,npleft,nprght,npplus,
		nelnsq,nother,ntotal,noobnd,nbevth;
    int		ev_min = MAXADU, xav = 0, yav = 0;
    short	min_adu = MAXADU, max_adu = 0;
    short	min_2ct = MAXADU, max_2ct = 0;
    short	xn = 512, xx = 0, yn = 512, yx = 0;

    int nacc, nnoto;
    unsigned char accmap[256], notomap[256], filter=0xff;
}

/*
 *  Initialize the histogram tables.  Each entry of the table needs
 *  to know which counter to increment, which histogram to increment
 *  and which extra pixels should be included in the summed pha,
 *  which only occurs for the L, Q and Other grades.
 */
static void
prep_hist()
{
    static
    const int extra[][4] = {  {4,4,4,4},
                              {0,4,4,4}, {2,4,4,4}, {6,4,4,4}, {8,4,4,4},
                              {0,2,4,4}, {0,6,4,4}, {6,8,4,4}, {8,2,4,4},
                              {0,2,6,8} };
    static
    const unsigned char	emask[] = { 0x00,
			    0x0a,0x12,0x48,0x50,
			    0x1a,0x4a,0x58,0x52,
			    0x5a };
    int		histo[8][MAXADU];
    unsigned char	sngle[] = { 0x00 }; 					/* 0 */
    unsigned char	splus[] = { 0x01,0x04,0x20,0x80,0x05,0x21,0x81,		/* 1 */
                                    0x24,0x84,0xa0,0x25,0x85,0xa4,0xa1,0xa5 };
    unsigned char	pvert[] = { 0x02,0x40,0x22,0x82,0x41,0x44,0x45,0xa2 };	/* 2 */
    unsigned char	pleft[] = { 0x08,0x0c,0x88,0x8c };			/* 3 */
    unsigned char	prght[] = { 0x10,0x30,0x11,0x31 };			/* 4 */
/*
 * unsigned char	phorz[] = { 0x08,0x10,0x0c,0x88,0x30,0x11,0x8c,0x31 };
 */
//unsigned char	pplus[] = { 0xc0,0x60,0x90,0x14,0x03,0x06,0x09,0x28,	/* 5 */
//			    0x86,0x83,0x15,0x0d,0x89,0x8d,0xc4,0xc1,
//			    0xc5,0x26,0xa6,0x23,0xa3,0x34,0xa8,0x35,
//			    0xa9,0x2c,0xa8,0xac,0x64,0x61,0x65 };

    unsigned char	pplus[] = { 0x03,0x06,0x09,0x28,0x60,0xc0,0x90,0x14,	/* 5 */
                                    0x83,0x26,0x89,0x2c,0x64,0xc1,0x91,0x34,
                                    0x23,0x86,0x0d,0xa8,0x61,0xc4,0xb0,0x15,
                                    0xa3,0xa6,0x8d,0xac,0x65,0xc5,0xb1,0x35 };
/*
 *  gned char	elish[] = { 0x12,0x32,0x50,0x51,0x48,0x4c,0x0a,0x8a };
 *  gned char	squar[] = { 0x16,0xd0,0x68,0x0b,0x36,0xd1,0x6c,0x8b };
 */
    unsigned char	elnsq[] = { 0x12,0x32,0x50,0x51,0x48,0x4c,0x0a,0x8a,	/* 6 */
                                    0x16,0xd0,0x68,0x0b,0x36,0xd1,0x6c,0x8b };

    unsigned char	other[256] /* = { all of the rest } */;			/* 7 */

    /* zero everything in sight */
    nsngle = nsplus = npvert = npleft = nprght = npplus =
        nelnsq = nother = ntotal = noobnd = nbevth = 0;
    bzero((char *)histo, sizeof(histo));

    /* make up the array of acceptable maps. */

    for (int i = 0; i < 256; i++) accmap[i]=0;
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
    for (int i = 0; i < sizeof(other); i++) {
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
}

/*
 *  Accumulate the num events in the tables
 */
static int
make_hist(int event,
          int split,
          struct data_str *ev,
          double rst,
          char *sty,
          int phlo,
          int phhi          
         )
{
    /*
     *  Get some gross event parameters
     */
    if (ev->data[4] < ev_min) ev_min = ev->data[4];
    if (ev->data[4] < event) { nbevth++; return 0; }
    /*
     *  Insert the reset clock correction
     */
    switch (*sty) {
      case '6':
        ev->data[7] -= ev->data[6]*rst;
        ev->data[4] -= ev->data[3]*rst;
        ev->data[1] -= ev->data[0]*rst;
      case '3':
        ev->data[8] -= ev->data[7]*rst;
        ev->data[2] -= ev->data[1]*rst;
      case '1':
        ev->data[5] -= ev->data[4]*rst;
    }
    /*
     *  Characterize event & accumulate most of pha
     */
    unsigned char	map = 0;
    short sum = 0;
    short phe[9];
                
    for (int j = 0; j < 9; j++) {
        const short phj = ev->data[j];
        phe[j] = phj;
        if (phj < split && j != 4) {
            phe[j]=0;
            continue;
        }
        switch (j) {
          case 0: map |= 0x01; //phe[0]=phj; 
            break;
          case 1: map |= 0x02; sum += phj; break;
          case 2: map |= 0x04; //phe[2]=phj; 
            break;
          case 3: map |= 0x08; sum += phj; break;
          case 4: 	       sum += phj; break;
          case 5: map |= 0x10; sum += phj; break;
          case 6: map |= 0x20; //phe[6]=phj;
            break;
          case 7: map |= 0x40; sum += phj; break;
          case 8: map |= 0x80; //phe[8]=phj;
            break;
        }
    }

    /* 
     *  grade is identified. check with filter to see whether
     *  to pass it on or not.
     *
     */

    for (int j = 0; j < nacc;j++) {
        if (map == accmap[j]) {
            goto accept_event;
        }
    }
    if (filter & 0x80) {
        for (int j = 0; j< nnoto; j++) {
            if (map == notomap[j]) goto nextevent;
        }
        goto accept_event;
    }
nextevent: 
    return 0;
accept_event:
    /*
     *  Finish pha with extra pixels of L, Q, and O events
     */
    look_up *ent = &table[map];
    const int *xtr = ent->extr;
    for (int j = 0; xtr[j] != 4 && j < 4; j++) sum += phe[xtr[j]];
    /*
     *  Accumulate statistics and various bounds
     */
    if (sum >= MAXADU) { noobnd++;  return 0; }
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

    return (sum >= phlo && sum < phhi) ? 1 : 0;
}

#if defined(MAIN)
/*
 *  Usage complaint message
 */
static void
usage()
{
	(void)fprintf(stderr, "Usage:  rv_gradefilt event split ");
        (void)fprintf(stderr, "-g glist... < evlist > evlist2\n\n");
        (void)fprintf(stderr, "\tevent   == event threshold\n");
        (void)fprintf(stderr, "\tsplit   == split threshold\n");
        (void)fprintf(stderr, "\tglist   == list of grades to pass\n");
        (void)fprintf(stderr, "\tevlist  == rv_style event list\n");
        (void)fprintf(stderr, "\tevlist2 == filtered rv_style event list\n");
	(void)fprintf(stderr, "\nThis version uses Dec92 exclusive grades:\n");
	(void)fprintf(stderr, "\n\tNew Grade Definition for Bright Mode\n");
	(void)fprintf(stderr, "\t\n");
	(void)fprintf(stderr, "\tgrade    name       description\n");
	(void)fprintf(stderr, "\t\n");
	(void)fprintf(stderr, "\t0        single             pure single\n");
	(void)fprintf(stderr, "\t1        single+            single + corner(s)\n");
	(void)fprintf(stderr, "\t2        vertical split     vertical split (+detouched corner(s))\n");
	(void)fprintf(stderr, "\t3        left split         left  split (+detouched corner(s))\n");
	(void)fprintf(stderr, "\t4        right split        right split (+detouched corner(s))\n");
	(void)fprintf(stderr, "\t5        single-sided+      single-sided split + a touched corner\n");
	(void)fprintf(stderr, "\t6        L+square           L or square (+a detouched corner)\n");
	(void)fprintf(stderr, "\t7        Others             all others\n");
}

int
main(int argc, char **argv)
{
	int	event, split, num, gr, tot = 0;
        int     phlo=0, phhi=MAXADU - 1;
	char	style[256];
	double	reset=0;

	if (argc<2) {
	  usage();
	  return(1);
	}

	--argc;argv++;	event=atoi(*argv);
	--argc;argv++;	split=atoi(*argv);

	while (--argc>0) {
	  argv++;
	  switch (argv[0][0]) {
	  case '-':
	    switch (argv[0][1]) {
	    case 'g':
              {
                  unsigned char   grades[] = { 0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80 };

                  if (filter == 0xff) filter=0x00;
                  while (--argc>0) {
                      argv++;
                      if (argv[0][0] == '-') {
                          argc++;  argv--;  goto end_of_grades;
                      }
                      gr=(int)atoi(*argv);
                      if ((gr >= 8) || (gr < 0)) { usage(); return (1);}
                      filter |= grades[gr];
                  }
              end_of_grades:;
              }
	      break;
	    case 'p':
	      if (argc<2) {
		usage();
		return(1);
	      }
	      --argc;argv++;	phlo=atoi(*argv);
	      --argc;argv++;	phhi=atoi(*argv);
	      break;
	    default:
	      break;
	    }
	    break;
	  default:
	    break;
	  }
	}

        /* ready for the data now. */
	prep_hist();
        const int EVENTS = 1024;
        data_str eventdata[EVENTS];
	while ((num = fread((void *)eventdata, sizeof(data_str), EVENTS, stdin)) > 0) {
            tot += num;
            for (data_str *ev = eventdata; ev != eventdata + num; ++ev) {
                if (make_hist(event, split, ev, reset, style, phlo, phhi)) {
                    fwrite(ev, sizeof(data_str), 1, stdout);
                }
            }
        }

	return(0);
}
#endif
