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
        HistogramTableBase::RESET_STYLES style = HistogramTableBase::TNONE;
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

        HistogramTableBase table(event, split, style, reset);
        table.setFilter(filter);

        if (phhi < 0) {
            phhi = HistogramTableBase::MAXADU;
        }

        const int EVENTS = 1024;
        data_str eventdata[EVENTS];
	while ((num = fread((void *)eventdata, sizeof(data_str), EVENTS, stdin)) > 0) {
            tot += num;
            for (data_str *ds = eventdata; ds != eventdata + num; ++ds) {
                lsst::rasmussen::Event ev(*ds);

                if (table.process_event(&ev)) {
                    if (table.sum >= phlo && table.sum < phhi) {
                        fwrite(&ev, sizeof(data_str), 1, stdout); // n.b. slices the Event back to data_str
                    }
                }
            }
        }

	return(0);
}
#endif
