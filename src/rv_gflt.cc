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
#include <limits>
#include <algorithm>
#include "lsst/rasmussen/tables.h"

/*
 *  Accumulate the num events in the tables
 */
bool
HistogramTableGflt::process_event(const data_str *ev
                                 )
{
    /*
     *  Get some gross event parameters
     */
    if (ev->data[4] < ev_min) ev_min = ev->data[4];
    if (ev->data[4] < _event) {
        nbevth++;
        grd = -1;                       // We don't know map yet. N.b. grd is in HistogramTableBase.
        return false;
    }

    short phe[9];
    std::copy(ev->data, ev->data + 9, phe);

    applyResetClockCorrection(phe);
    /*
     *  Characterize event & accumulate most of pha
     */
    unsigned char map = 0;

    sum = 0;                            // in HistogramTableBase
    for (int j = 0; j < 9; j++) {
        const short phj = phe[j];

        if (phj < _split && j != 4) {
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
     *  grade is identified. check with _filter to see whether
     *  to pass it on or not.
     */

    bool accept = false;
    for (int j = 0; j < nacc;j++) {
        if (map == accmap[j]) {
            accept = true;
            break;
        }
    }
    if (!accept && (_filter & 0x80)) {
        for (int j = 0; j < nnoto; j++) {
            if (map == notomap[j]) {
                grd = setGrdFromType(map); // n.b. grd is in HistogramTableBase
                return false;
            }
        }
        accept = true;
    }

    if (!accept) {
        grd = setGrdFromType(map);      // n.b. grd is in HistogramTableBase
        return false;
    }

    return finishEventProcessing(ev, phe, map);
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
        int     phlo=0, phhi=-1;
        HistogramTableGflt::RESET_STYLES style = HistogramTableGflt::TNONE;
	double	reset=0;

	if (argc<2) {
	  usage();
	  return(1);
	}

	--argc;argv++;	event=atoi(*argv);
	--argc;argv++;	split=atoi(*argv);

        int filter = 0x00;              // mask of grades we care about

	while (--argc>0) {
	  argv++;
	  switch (argv[0][0]) {
	  case '-':
	    switch (argv[0][1]) {
	    case 'g':
              {
                  unsigned char   grades[] = { 0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80 };

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

        HistogramTableGflt table(filter, event, split, style, reset);

        if (phhi < 0) {
            phhi = HistogramTableGflt::MAXADU;
        }

        const int EVENTS = 1024;
        data_str eventdata[EVENTS];
	while ((num = fread((void *)eventdata, sizeof(data_str), EVENTS, stdin)) > 0) {
            tot += num;
            for (data_str *ev = eventdata; ev != eventdata + num; ++ev) {
                if (table.process_event(ev)) {
                    if (table.sum >= phlo && table.sum < phhi) {
                        fwrite(ev, sizeof(data_str), 1, stdout);
                    }
                }
            }
        }

	return(0);
}
#endif
