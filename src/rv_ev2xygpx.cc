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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits>
#include <algorithm>
#include "lsst/rasmussen/tables.h"

/*
 *  Accumulate the num events in the tables
 */
bool
HistogramTableXygpx::process_event(lsst::rasmussen::Event *ev
                                  )
{
    /*
     *  Get some gross event parameters
     */
    if (ev->data[4] < ev_min) ev_min = ev->data[4];
    if (ev->data[4] < _event) { nbevth++; return false; }

    short phe[9];
    std::copy(ev->data, ev->data + 9, phe);

    applyResetClockCorrection(phe);
    /*
     *  Characterize event & accumulate most of pha
     */
    unsigned char map = 0;

    p9 = 0;                             // in HistogramTableBase
    sum = 0;                            // in HistogramTableBase
    for (int j = 0; j < 9; j++) {
        const short phj = phe[j];

        switch (_do_what) {
          case p_9:
            p9 += phj;
            break;
          case p_1357:
            if (j == 1 || j == 3 || j == 5 || j == 7 || j == 4) p9 += phj;
            break;
          case p_17:
            if (j == 1 || j == 7 || j == 4) p9 += phj;
            break;
          case p_35:
            if (j == 3 || j == 5 || j == 4) p9 += phj;
            break;
          case p_list:
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
    return finishEventProcessing(ev, phe, map);
}

#if defined(MAIN)
/*
 *  Usage complaint message
 */
static void
usage()
{
	(void)fprintf(stderr, "Usage:  rv_ev2pcf [--ev2pcf] event split sfile ");
	(void)fprintf(stderr, "[reset [style]] < evlist > pcfile\n\n");
	(void)fprintf(stderr, "\tevent  == event threshold\n");
	(void)fprintf(stderr, "\tsplit  == split threshold\n");
	(void)fprintf(stderr, "\tsfile  == source/expt data file\n");
	(void)fprintf(stderr, "\treset  == reset correction factor\n");
	(void)fprintf(stderr, "\tstyle  == of correction: 1 | 3 | 6\n");
	(void)fprintf(stderr, "\tevlist == rv_style event list\n");
	(void)fprintf(stderr, "\tpcfile == qdp_histogram/pcf data\n");

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
	(void)fprintf(stderr, "\n");
	(void)fprintf(stderr, "With no arguments, a CLASSIFY format\n");
	(void)fprintf(stderr, "lookup grade table is dumped to stderr.\n");
}

/*
 * Write out the x,y,grade,ph for each acceptable event.
 */
int
main(int argc, char **argv)
{
    const bool ev2pcf = (argc > 1 && strcmp(argv[1], "--ev2pcf") == 0); // run the ev2pcf code
    if (ev2pcf) {
        --argc; ++argv;
    }

	int	event, split, num, tot = 0;
        HistogramTableXygpx::RESET_STYLES style = HistogramTableXygpx::T1;

	const char *sfile = "unknown";
	double	reset = 0;
	char    *calc;

	if (argc == 1) {	/* for diagnostic purposes */
            HistogramTableXygpx table;
            table.dump_table();
            return 0;
	}
	if (argc < 3 || argc > 6) { usage(); return 1; }

        HistogramTableXygpx::calctype do_what = ev2pcf ? HistogramTableXygpx::p_list : HistogramTableXygpx::p_9;
	if (--argc > 0) {
	  event = atoi(*++argv);
	  if (--argc > 0) {
	    split = atoi(*++argv);
	    if (--argc > 0) {
	      calc=*++argv;
	      if (strcmp(calc,"p9")==0) {
                  do_what=HistogramTableXygpx::p_9;
	      } else if (strcmp(calc,"p17")==0) {
		  do_what=HistogramTableXygpx::p_17;
              } else if (strcmp(calc,"p35")==0) {
                  do_what=HistogramTableXygpx::p_35;
              } else if (strcmp(calc,"p1357")==0) {
                  do_what=HistogramTableXygpx::p_1357;
              } else if (strcmp(calc,"plist")==0) {
                  do_what=HistogramTableXygpx::p_list;
              } else {
                  fprintf(stderr,"don't recognize this arg: %s\n exiting..",calc);
                  exit(1);
	      }
	      if (--argc) {
		sfile = *++argv;
		if (--argc) {
		  reset = atof(*++argv);
		  if (--argc) {
                      const char *styleStr = *++argv;
                                            
                      switch (*styleStr) {
                        case '1': style = HistogramTableXygpx::T1; break;
                        case '3': style = HistogramTableXygpx::T3; break;
                        case '6': style = HistogramTableXygpx::T6; break;
                        default:
                          fprintf(stderr,"Invalid reset style: %s\n exiting..", styleStr);
                          return 1;
                      }
		  }
		}
	      }
	    }
	  }
	}

        HistogramTableXygpx table(do_what, event, split, style, reset);
        const int EVENTS = 1024;
        data_str eventdata[EVENTS];

	while ((num = fread((void *)eventdata, sizeof(data_str), EVENTS, stdin)) > 0) {
            tot += num;
            for (data_str *ds = eventdata; ds != eventdata + num; ++ds) {
                lsst::rasmussen::Event ev(*ds);
                if (table.process_event(&ev)) {
                    if (ev2pcf) {
                        continue;
                    }
                    if (do_what == HistogramTableXygpx::p_list) {
                        fprintf(stdout,"%d %d %d %d p:", ev.x, ev.y, ev.grade, table.sum);
                        for (int i=0;i<9;i++) {
                            fprintf(stdout," %f",ev.data[i]);
                        }
                        fprintf(stdout,"\n");
                    } else {
                        fprintf(stdout,"%d %d %d %d %g %d\n",ev.x,ev.y, ev.grade, table.sum,
                                ev.data[4], table.p9);
                    }
                    
                }
            }
	}

        if (ev2pcf) {
            table.dump_head(stdout, sfile, tot);
            table.dump_hist(stdout, sfile);
        }

	return 0;
}
#endif
