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

namespace {
class HistogramTable : public HistogramTableBase {
public:
    enum calctype { p_9,
                    p_17,
                    p_35,
                    p_1357,
                    p_list,             // for the "total"
    };

    HistogramTable(calctype do_what=p_list) : HistogramTableBase(), _do_what(do_what) {}
    virtual int process_event(const data_str *ev, int event, int split, RESET_STYLES sty, double rst);

    // Value set by process_event
    int p9;

private:
    const calctype _do_what;
};

/*
 *  Accumulate the num events in the tables
 */
int
HistogramTable::process_event(
        const data_str *ev,
        int event,
        int split,
        RESET_STYLES sty,
        double rst
                             )
{
    /*
     *  Get some gross event parameters
     */
    if (ev->data[4] < ev_min) ev_min = ev->data[4];
    if (ev->data[4] < event) { nbevth++; return 0; }

    short phe[9];
    std::copy(ev->data, ev->data + 9, phe);
    /*
     *  Insert the reset clock correction
     */
    switch (sty) {
      case T6:
        phe[7] -= phe[6]*rst;
        phe[4] -= phe[3]*rst;
        phe[1] -= phe[0]*rst; // fall through
      case T3:
        phe[8] -= phe[7]*rst;
        phe[2] -= phe[1]*rst; // fall through
      case T1:
        phe[5] -= phe[4]*rst;
        break;
      case TNONE:
        break;
    }
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

        if (phj < split && j != 4) {
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
}


/*
 *  Dump the basic calibration file header
 */
const int NAMLEN = 512;

static void
dump_head(const HistogramTable &table,
          const char *sfile, char *efile,
          int event, int split, int total)          
{
	FILE	*fp;
	char	line[NAMLEN], *c;

	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  QDP Basic Calibration File\n");
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Working_dir  = %s\n", getcwd((char *)0, NAMLEN));
	(void)fprintf(stdout, "!  Source_file  = %s\n", sfile);
	(void)fprintf(stdout, "!  Event_thresh = %d\n", event);
	(void)fprintf(stdout, "!  Split_thresh = %d\n", split);
	//(void)fprintf(stdout, "!  Total_frames = %d\n", table.cnt);
	(void)fprintf(stdout, "!  Total_events = %d\n", table.ntotal);
	(void)fprintf(stdout, "!  Total_pixels = %d\n", (table.xx - table.xn)*(table.yx - table.yn));
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Events_below = %d\n", table.nbevth);
	(void)fprintf(stdout, "!  Events_above = %d\n", table.noobnd);
	(void)fprintf(stdout, "!  Events_input = %d\n", total);
	(void)fprintf(stdout, "!  PH4_minimum  = %d\n", table.ev_min);
	(void)fprintf(stdout, "!  PHS_minimum  = %d\n", table.min_adu);
	(void)fprintf(stdout, "!  PHS_Maximum  = %d\n", table.max_adu);
	(void)fprintf(stdout, "!  X_Minimum    = %d\n", table.xn);
	(void)fprintf(stdout, "!  X_Average    = %d\n", table.xav/table.ntotal);
	(void)fprintf(stdout, "!  X_Maximum    = %d\n", table.xx);
	(void)fprintf(stdout, "!  Y_Minimum    = %d\n", table.yn);
	(void)fprintf(stdout, "!  Y_Average    = %d\n", table.yav/table.ntotal);
	(void)fprintf(stdout, "!  Y_Maximum    = %d\n", table.yx);

	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Exclusive grades -- corrected L+Q.\n");
	(void)fprintf(stdout, "!\n");

        if (strcmp(sfile, "unknown") == 0) {
            return;
        }
        
	if ((fp = fopen(sfile, "r")) == NULL) return;
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Experimental parameters\n");
	(void)fprintf(stdout, "!\n");
	while (fgets(line, NAMLEN-1, fp)) {
		if (line[0] == '#') continue;
		(void)fprintf(stdout, "!  %s", line);
		if (!strncmp(line, "evlist  = ", 10)) {
			(void)strcpy(efile, &line[10]);
			for (c = efile; *c; c++)
				if (*c == '\t') *c = ' '; 
		}
	}
	(void)fclose(fp);
}

#if defined(MAIN)
/*
 *  Usage complaint message
 */
static void
usage()
{
	(void)fprintf(stderr, "Usage:  rv_ev2pcf event split sfile ");
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
        HistogramTable::RESET_STYLES style = HistogramTable::T1;

	const char *sfile = "unknown";
	double	reset = 0;
	char    *calc;

	if (argc == 1) {	/* for diagnostic purposes */
            HistogramTable table;
            table.dump_table();
            return 0;
	}
	if (argc < 3 || argc > 6) { usage(); return 1; }

        HistogramTable::calctype do_what = ev2pcf ? HistogramTable::p_list : HistogramTable::p_9;
	if (--argc > 0) {
	  event = atoi(*++argv);
	  if (--argc > 0) {
	    split = atoi(*++argv);
	    if (--argc > 0) {
	      calc=*++argv;
	      if (strcmp(calc,"p9")==0) {
                  do_what=HistogramTable::p_9;
	      } else if (strcmp(calc,"p17")==0) {
		  do_what=HistogramTable::p_17;
              } else if (strcmp(calc,"p35")==0) {
                  do_what=HistogramTable::p_35;
              } else if (strcmp(calc,"p1357")==0) {
                  do_what=HistogramTable::p_1357;
              } else if (strcmp(calc,"plist")==0) {
                  do_what=HistogramTable::p_list;
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
                        case '1': style = HistogramTable::T1; break;
                        case '3': style = HistogramTable::T3; break;
                        case '6': style = HistogramTable::T6; break;
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

        HistogramTable table(do_what);
        const int EVENTS = 1024;
        data_str eventdata[EVENTS];

	while ((num = fread((void *)eventdata, sizeof(data_str), EVENTS, stdin)) > 0) {
            tot += num;
            for (data_str *ev = eventdata; ev != eventdata + num; ++ev) {
                if (table.process_event(ev, event, split, style, reset)) {
                    if (ev2pcf) {
                        continue;
                    }
                    if (do_what == HistogramTable::p_list) {
                        fprintf(stdout,"%d %d %d %d p:", ev->x, ev->y, table.grd, table.sum);
                        for (int i=0;i<9;i++) {
                            fprintf(stdout," %f",ev->data[i]);
                        }
                        fprintf(stdout,"\n");
                    } else {
                        fprintf(stdout,"%d %d %d %d %g %d\n",ev->x,ev->y, table.grd, table.sum,
                                ev->data[4], table.p9);
                    }
                    
                }
            }
	}

        if (ev2pcf) {
            char efile[NAMLEN];
            dump_head(table, sfile, efile, event, split, tot);
            table.dump_hist(event, split, sfile, "");
        }

	return 0;
}
#endif
