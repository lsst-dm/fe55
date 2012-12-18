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
#include <unistd.h>
#include "lsst/rasmussen/rv.h"

#define EVENTS 		1024
#define MAXADU 		4096
#define EXTADU		8
#define NAMLEN		512

struct data_str eventdata[EVENTS];
struct look_up { int *type, *extr, *hist; } table[256];

int		nsngle,nsplus,npvert,npleft,nprght,npplus,
		nelnsq,nother,ntotal,noobnd,nbevth;
int		extra[][4] = {  {4,4,4,4},
				{0,4,4,4}, {2,4,4,4}, {6,4,4,4}, {8,4,4,4},
				{0,2,4,4}, {0,6,4,4}, {6,8,4,4}, {8,2,4,4},
				{0,2,6,8} };
unsigned char	emask[] = { 0x00,
			    0x0a,0x12,0x48,0x50,
			    0x1a,0x4a,0x58,0x52,
			    0x5a };
int		histo[8][MAXADU], oldfrnum = 0,
		cnt = 0, ev_min = MAXADU, xav = 0, yav = 0;
short		min_adu = MAXADU, max_adu = 0;
short		min_2ct = MAXADU, max_2ct = 0;
short		xn = 512, xx = 0, yn = 512, yx = 0;

unsigned char	sngle[] = { 0x00 }; 					/* 0 */
unsigned char	splus[] = { 0x01,0x04,0x20,0x80,0x05,0x21,0x81,		/* 1 */
			    0x24,0x84,0xa0,0x25,0x85,0xa4,0xa1,0xa5 };
unsigned char	pvert[] = { 0x02,0x40,0x22,0x82,0x41,0x44,0x45,0xa2 };	/* 2 */
unsigned char	pleft[] = { 0x08,0x0c,0x88,0x8c };			/* 3 */
unsigned char	prght[] = { 0x10,0x30,0x11,0x31 };			/* 4 */
/*
 *  gned char	phorz[] = { 0x08,0x10,0x0c,0x88,0x30,0x11,0x8c,0x31 };
 */
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
char		efile[NAMLEN];

static void	usage(), prep_hist(), make_classification(),
		dump_head(), dump_hist(), dump_table();

typedef enum calctype { p_9,
			p_17,
			p_35,
			p_1357,
			p_list } calctype; // for the "total"
calctype do_what=p_9;

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
 *  For diagnostic purposes, dump the grade table in CLASSIFY format.
 */
static void
dump_table()
{
	register int		i, j;

	for (i = 0; i < 256; i++) {

		if ( table[i].type == &nsngle )  j = 0;	else
		if ( table[i].type == &nsplus )  j = 1;	else
		if ( table[i].type == &npvert )  j = 2;	else
		if ( table[i].type == &npleft )  j = 3;	else
		if ( table[i].type == &nprght )  j = 4;	else
		if ( table[i].type == &npplus )  j = 5;	else
		if ( table[i].type == &nelnsq )  j = 6;	else
		if ( table[i].type == &nother )  j = 7;	else
						 j = 9;

		(void)fprintf(stderr, "%d,", j);
		if (i%16 == 15) (void)fprintf(stderr, "\n");
	}
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
	register int		i;
	register struct look_up	*t;
	unsigned char		b;
	int			j;

	/* zero everything in sight */
	nsngle = nsplus = npvert = npleft = nprght = npplus =
		 nelnsq = nother = ntotal = noobnd = nbevth = 0;
	bzero((char *)histo, sizeof(histo));

	/* load the sngle events into table GRADE 0 */
	for (i = 0; i < sizeof(sngle); i++) {
		t = table + sngle[i];
		t->type = &nsngle;
		t->hist = histo[0];
		t->extr = extra[0];	
	}

	/* load the splus events into table GRADE 1 */
	for (i = 0; i < sizeof(splus); i++) {
		t = table + splus[i];
		t->type = &nsplus;
		t->hist = histo[1];
		t->extr = extra[0];	
	}

	/* load the pvert events into table GRADE 2 */
	for (i = 0; i < sizeof(pvert); i++) {
		t = table + pvert[i];
		t->type = &npvert;
		t->hist = histo[2];
		t->extr = extra[0];	
	}

	/* load the pleft events into table GRADE 3 */
	for (i = 0; i < sizeof(pleft); i++) {
		t = table + pleft[i];
		t->type = &npleft;
		t->hist = histo[3];
		t->extr = extra[0];	
	}

	/* load the prght events into table GRADE 4 */
	for (i = 0; i < sizeof(prght); i++) {
		t = table + prght[i];
		t->type = &nprght;
		t->hist = histo[4];
		t->extr = extra[0];	
	}

	/* load the pplus events into table GRADE 5 */
	for (i = 0; i < sizeof(pplus); i++) {
		t = table + pplus[i];
		t->type = &npplus;
		t->hist = histo[5];
		t->extr = extra[0];	
	}

	/* load the elnsq events into table GRADE 6 */
	for (i = 0; i < sizeof(elnsq); i++) {
		t = table + elnsq[i];
		t->type = &nelnsq;
		t->hist = histo[6];
		t->extr = extra[0];	
		for (b = 0x5a & elnsq[i], j = 0; j < sizeof(emask); j++)
			if (b == emask[j]) {
				t->extr = extra[j];
				break;
			}
	}

	/* load the other events into table GRADE 7 */
	for (i = 0; i < sizeof(other); i++) {
		t = table + i;
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
static void
make_classification(int event,
                    int split,
                    int num,
                    struct data_str *ev,
                    double rst,
                    char *sty
        )
{
  // in this routine we write out the x,y,grade,ph for each event.
	register unsigned char	map;
	register int		j;
	short			phj, sum, phe[9], hsum;
	struct look_up		*ent;
	int			*xtr;

	for ( ; num--; ev++) {
		/*
		 *  Get some gross event parameters
		 */
		if (oldfrnum != ev->framenum) {
			oldfrnum = ev->framenum;
			cnt++;
		}
		if (ev->data[4] < ev_min) ev_min = ev->data[4];
		if (ev->data[4] < event) { nbevth++; continue; }
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
			default:	break;		/* no correction */
		}
		/*
		 *  Characterize event & accumulate most of pha
		 */
		int p4 = ev->data[4];
		int p9 = 0;
		for (j = 0, map = 0, sum = 0; j < 9; j++) {
			phe[j] = phj = ev->data[j];

			if (do_what==p_9) {
			  p9 += phj;
			} 
			if (do_what==p_1357) {
			  if ((j-1)*(j-3)*(j-5)*(j-7)*(j-4)==0) p9+=phj;
			}
			if (do_what==p_17) {
			  if ((j-1)*(j-7)*(j-4)==0) p9+=phj;
			}
			if (do_what==p_35) {
			  if ((j-3)*(j-5)*(j-4)==0) p9+=phj;
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
				case 4: 	     sum += phj; break;
				case 5: map |= 0x10; sum += phj; break;
				case 6: map |= 0x20;           ; break;
				case 7: map |= 0x40; sum += phj; break;
				case 8: map |= 0x80;           ; break;
			}
		}
		/*
		 *  Finish pha with extra pixels of L, Q, and O events
		 */
		ent = &table[map];
		xtr = ent->extr;
		for (j = 0; xtr[j] != 4 && j < 4; j++) sum += phe[xtr[j]];
		/*
		 *  Accumulate statistics and various bounds
		 */
		if (sum >= MAXADU) { noobnd++;  continue; }
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
		hsum = ent->hist[sum] += 1;
		if (hsum > 2) {
			if (sum > max_2ct) max_2ct = sum;
			if (sum < min_2ct) min_2ct = sum;
		}
		{
		  int grd;
		  if (ent->type == &nsngle) { grd=0; } 
		  else { if (ent->type == &nsplus) { grd=1; }
		    else { if (ent->type == &npvert) {	grd=2; }
		      else { if (ent->type == &npleft) { grd=3; }
			else { if (ent->type == &nprght) { grd=4; } 
			  else { if (ent->type == &npplus) { grd=5; } 
			    else { if (ent->type == &nelnsq) { grd=6; } 
			      else { if (ent->type == &nother) { grd=7; }
				else { grd=-1; } } } } } } } }
		  if (do_what == p_list) {
		    fprintf(stdout,"%d %d %d %d p:",ev->x,ev->y,grd,sum);
		    {
		      int i;
		      for (i=0;i<9;i++) 
			fprintf(stdout," %g", ev->data[i]);
		      fprintf(stdout,"\n");
		    }
		  } else {
		    fprintf(stdout,"%d %d %d %d %d %d\n",ev->x,ev->y,grd,sum,p4,p9);
		  }
		}
	}
}

/*
 *  Dump the basic calibration file header
 */
static void
dump_head(char *sfile, int event, int split, int total)
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
	(void)fprintf(stdout, "!  Total_frames = %d\n", cnt);
	(void)fprintf(stdout, "!  Total_events = %d\n", ntotal);
	(void)fprintf(stdout, "!  Total_pixels = %d\n", (xx - xn)*(yx - yn));
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Events_below = %d\n", nbevth);
	(void)fprintf(stdout, "!  Events_above = %d\n", noobnd);
	(void)fprintf(stdout, "!  Events_input = %d\n", total);
	(void)fprintf(stdout, "!  PH4_minimum  = %d\n", ev_min);
	(void)fprintf(stdout, "!  PHS_minimum  = %d\n", min_adu);
	(void)fprintf(stdout, "!  PHS_Maximum  = %d\n", max_adu);
	(void)fprintf(stdout, "!  X_Minimum    = %d\n", xn);
	(void)fprintf(stdout, "!  X_Average    = %d\n", xav/ntotal);
	(void)fprintf(stdout, "!  X_Maximum    = %d\n", xx);
	(void)fprintf(stdout, "!  Y_Minimum    = %d\n", yn);
	(void)fprintf(stdout, "!  Y_Average    = %d\n", yav/ntotal);
	(void)fprintf(stdout, "!  Y_Maximum    = %d\n", yx);

	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Exclusive grades -- corrected L+Q.\n");
	(void)fprintf(stdout, "!\n");

	if ((fp = fopen(sfile, "r")) == NULL) return;
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Experimental parameters\n");
	(void)fprintf(stdout, "!\n");
	while (fgets(line, NAMLEN-1, fp)) {
		if (line[0] == '#') continue;
		fprintf(stdout, "!  %s", line);
		if (!strncmp(line, "evlist  = ", 10)) {
		  strcpy(efile, &line[10]);
		  for (c = efile; *c; c++)
		    if (*c == '\t') *c = ' '; 
		}
	}
	(void)fclose(fp);
}

/*
 *  Dump the QDP header histogram table
 */
static void
dump_hist(int event, int split, char *sfile)
{
	register int	i;

	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  QDP Header follows\n");
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "lab top Event = %d Split = %d Source = %s\n",
		event, split, sfile);
	if (efile[0]) (void)fprintf(stdout, "lab file %s", efile);
	(void)fprintf(stdout, "lab g1 Pulse Height (ADU)\n");
	(void)fprintf(stdout, "lab rot\n");
	(void)fprintf(stdout, "lab g2 N(S)\nlab g3 N(S+)\n");
	(void)fprintf(stdout, "lab g4 N(Pv)\nlab g5 N(Pl)\n");
	(void)fprintf(stdout, "lab g6 N(Pr)\nlab g7 N(P+)\n");
	(void)fprintf(stdout, "lab g8 N(L+Q)\nlab g9 N(O)\n");
	(void)fprintf(stdout, "csize 0.75\n");

	min_2ct = (  min_2ct <        EXTADU) ?      0 : min_2ct - EXTADU;
	max_2ct = (++max_2ct > MAXADU-EXTADU) ? MAXADU : max_2ct + EXTADU;
	(void)fprintf(stdout, "res x %d %d\n", min_2ct, max_2ct);
	(void)fprintf(stdout, "res y2 1\nres y3 1\n");
	(void)fprintf(stdout, "res y4 1\nres y5 1\n");
	(void)fprintf(stdout, "res y6 1\nres y7 1\n");
	(void)fprintf(stdout, "res y8 1\nres y9 1\n");

	(void)fprintf(stdout, "error y sq 2 3 4 5 6 7 8 9\n");
	(void)fprintf(stdout, "log y on\n");
	(void)fprintf(stdout, "plot vert\n");
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  Histogram data follows\n");
	(void)fprintf(stdout, "!\n");
	(void)fprintf(stdout, "!  PHA\tS\tS+\tPv\tPl\tPr\tP+\tL+Q\tO\n");
	(void)fprintf(stdout, "!  TOT\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		nsngle,nsplus,npvert,npleft,nprght,npplus,nelnsq,nother);
	(void)fprintf(stdout, "!\n");

	min_adu = (  min_adu <        EXTADU) ?      0 : min_adu - EXTADU;
	max_adu = (++max_adu > MAXADU-EXTADU) ? MAXADU : max_adu + EXTADU;
	for (i = min_adu; i < max_adu; i++)
		(void)fprintf(stdout,
			"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i,
			histo[0][i], histo[1][i], histo[2][i], histo[3][i],
			histo[4][i], histo[5][i], histo[6][i], histo[7][i]);
}


int
main(int argc, char **argv)
{
	int	event, split, num, tot = 0;
	char	*sfile, def_style = '1', *style = &def_style;
	double	reset = 0;
	char    *calc;

	if (argc == 1) {	/* for diagnostic purposes */
		prep_hist();
		dump_table();
		return(0);
	}
	if (argc < 3 || argc > 6) { usage(); return(1); }
	if (--argc) {
	  event = atoi(*++argv);
	  if (--argc) {
	    split = atoi(*++argv);
	    if (--argc) {
	      calc=*++argv;
	      if (strcmp(calc,"p9")==0) {
		do_what=p_9;
	      } else {
		if (strcmp(calc,"p17")==0) {
		  do_what=p_17;
		} else {
		  if (strcmp(calc,"p35")==0) {
		    do_what=p_35;
		  } else {
		    if (strcmp(calc,"p1357")==0) {
		      do_what=p_1357;
		    } else {
		      if (strcmp(calc,"plist")==0) {
			do_what=p_list;
		      } else {
			fprintf(stderr,"don't recognize this arg: %s\n exiting..",calc);
			exit(1);
		      }
		    }
		  }
		}
	      }
	      if (--argc) {
		sfile = *++argv;
		if (--argc) {
		  reset = atof(*++argv);
		  if (--argc) {
		    style = *++argv;
		  }
		}
	      }
	    }
	  }
	}

	prep_hist();
	int nread = 1000;		/* number of events to read at a time */
	while ((num = fread((void *)eventdata, sizeof(struct data_str), nread, stdin)) > 0) {
		tot += num;
		make_classification(event, split, num, eventdata, reset, style);
	}
	//	dump_head(sfile, event, split, tot);
	//	dump_hist(event, split, sfile);
	return(0);
}
