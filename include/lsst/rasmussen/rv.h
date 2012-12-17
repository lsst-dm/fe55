/* rv.h		24 July 89	Roland Vanderspek

   Include file for RV tools software

   Modified 15 April 1993.  Adjusted to simplify code for data analysis,
     removing the emphasis on DE emulation.  Dark frames disappear, along
     with DARK_FB_BITS and virtual frame buffers and other such stuff. 
     New data structure is gone, and perhaps the preamble will go, too.
*/

#define NUM_UNDERCLOCKS		4
#define NUM_OVERCLOCKS		18
#define NUM_POVERCLOCKS		20

#define NUM_FM_BIASROWS		8	/* # of FAST rows which use same bias */
#define NUM_TOSUM_FM_BIASROWS	1	/* # of FAST rows summed for bias */

#define SYNCMASK		0x8000
#define CHIPMASK		0x6000
#define CHIPSHFT		13

#define param_file		"FINDEVENTS.PARAM"

struct data_str
  {
  char mode;
  int framenum;
  short chipnum,x,y,data[9];
  };
short datastr_size = sizeof(struct data_str);

short pream_size = 0;


