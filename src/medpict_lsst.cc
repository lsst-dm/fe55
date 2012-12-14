#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include "/usr/local/heasoft/heacore/cfitsio/fitsio.h"
#include "/usr/local/heasoft/heacore/cfitsio/fitsio2.h"

#include "rv.h"

#define FNMAX     2000
#define OCHISTMAX 300
#define OCHISTOS  150
#define OCMAX     8
#define HDR_MAXRCDS 8

enum format {e2v_ccd250,lsst_sta_studycontract,bnl_e2v_studycontract,lbox,astd,berlin,hirefs,raw};

int comp_int(const void *a,const void *b);
void median(int x[],int n,int *xmed);
void printerror( int status);

typedef struct ev_s {
  struct data_str event;
  struct ev_s *prev_ev;
} ev_stack;

char filename[FNMAX][1024]; 
long filepos[FNMAX]; /* filepos[fi] is 0L when the file is open, IF it is
		        to be closed later. [filepos shouldn't be used
			unless trying to emulate `unlimited' number of 
			descriptors]. If file must be closed, filepos[fi]
			should be used to record the last position 
			before closing.*/

main(int argc,char *argv[])
{
  char cmd[204800];
  char biasout[1024],header[2880*HDR_MAXRCDS+8],
       destination_directory[1024],eventfilename[FNMAX][1024],*occLUTAB,
       histfilename[FNMAX][1024],tempstr[1024],input_biasfile[1024];
  ev_stack *evstack[FNMAX],*current_evstack_pointer;
  enum format ocspec;
  //  FILE *fp[FNMAX],*fib;
  FILE *fhist[FNMAX],*fevlist[FNMAX];
  int  status=0;
  fitsfile *ffp[FNMAX],*ffib=NULL,*ffout=NULL;

  int fni=0,fi,i,ti,pi,ncount,nx,ny,npix,fill,burstmode,nev,tmpint,
       *hist[FNMAX],eventsearch,x,y,xi,yi,evthresh,histmin,histmax,histmode,
       nocpix[OCMAX],ncor[OCMAX],ocsample[OCMAX][2048],ocsample_y[2048],
       OCcorrect[OCMAX][2048],nysample,oc,OCcorrection,tx,noc;

  int ocsLU[2048],occLU[2048];

  float PIX[FNMAX],med,weight[FNMAX+1],sx,sxx,mean,sdev,plim[2],
        arg,ocval[FNMAX][OCMAX];
  int pix[FNMAX],tmp,*medianbias,*cpix,ib,use_bias,
        oc_int[FNMAX][OCMAX];
  int *dp[FNMAX];
  long memorysize;
  struct data_str event;

  int REB=1;

  OCcorrection=0;  burstmode=0;               
  eventsearch= 0;  histmode=0;                histmin=-4000;
  evthresh=    20;destination_directory[0]=0;histmax=4000;
  biasout[0]=  0;  input_biasfile[0]=0;       medianbias=0;
  occLUTAB=    0;
  ocspec=bnl_e2v_studycontract;
  ocspec=lsst_sta_studycontract;

  // write entire cmdline into char cmd[]
  { 
    int i;
    sprintf(cmd,"%s",argv[0]);
    for (i=1;i<argc;i++) {
      sprintf(cmd,"%s %s",cmd,argv[i]);
    }
  }

  while (--argc > 0)  {
    argv++;
    if (argv[0][0]=='-') {
      switch (argv[0][1]) {
      case 'o':
	--argc;	argv++;
	sprintf(biasout,"%s",argv[0]);
	break;
      case 't':
	--argc; argv++;
	evthresh=(int)atof(argv[0]);
	fprintf(stderr,"setting event threshold to %d.\n",evthresh);
	break;
      case 'h':
	histmode=1;
	break;
      case 'R': // rebin specification
	--argc;argv++;
	REB=atoi(argv[0]);
	break;
      case 'f': /* format specification */
	--argc;argv++;
	switch (argv[0][0]) {
	case 'c':
	  ocspec=e2v_ccd250;
	  break;
	case 's':
	  ocspec=lsst_sta_studycontract;
	  break;
	case 'e':
	  ocspec=bnl_e2v_studycontract;
	  break;
	case 'l':
	  ocspec=lbox;
	  break;
	case 'h':
	  ocspec=hirefs;
	  break;
	case 'a':
	  ocspec=astd;
	  break;
	case 'b':
	  ocspec=berlin;
	  break;
	default:
	  usage("can't parse format specification.");
	}
	break;
      case 'b':
	burstmode=1;
	break;
      case 'e':
	eventsearch=1;
	break;
      case 'c':
	OCcorrection=1;
	break;
      case 'B':
	--argc; argv++;
	sprintf(input_biasfile,"%s",argv[0]);
	//	if ((fib=fopen(input_biasfile,"r"))==NULL) 
	//	  usage("can't open input bias file.");
	if (fits_open_file(&ffib,input_biasfile,READONLY,&status)) {
	  printerror(status);
	}
	break;
      case 'd':
	--argc; argv++;
	sprintf(destination_directory,"%s",argv[0]);
	if (strrchr(destination_directory,(int)'/')-destination_directory
	    != strlen(destination_directory)-1) {
	  sprintf(tempstr,"%s/",destination_directory);
	  sprintf(destination_directory,"%s",tempstr);
	}
	break;
      default:
	usage("unknown option.");
	break;
      }
    } else {
      /* it's a filename */
      //      fprintf(stderr,"got filename %s\n",argv[0]);
      sprintf(filename[fni],"%s",argv[0]);
      // old way: using FILE descriptors
      //      if ((fp[fni]=fopen(filename[fni],"r"))==NULL) usage("can't open file.");
      //      /* close the file again to keep options open */
      //      dp[fni]=0;
      //      fclose(fp[fni]);
      //      fp[fni]=NULL;
      //      filepos[fni]=0L;
      // 
      // new way: using fitsfiles
      if (fits_open_file(&ffp[fni],filename[fni],READONLY,&status))
	printerror(status);
      // close file again to keep options open
      dp[fni]=NULL;
      if (fits_close_file(ffp[fni],&status))
	printerror(status);
      ffp[fni]=NULL;
      filepos[fni]=1L;
      evstack[fni]=NULL;
      fni++;
    }
  }

  if (fni<FNMAX-5) { /* arbitrary limit - often csh:limit shows 64 descriptors */
    /* open the input files and leave them open */
    for (fi=0;fi<fni;fi++){
      if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	printerror(status);
      //      if ((fp[fi]=fopen(filename[fi],"r"))==NULL)
      //	usage("can't open file.");
      filepos[fi]=1L;
    }
  }

  if ( ! biasout[0] && ! eventsearch && ! histmode ) 
    usage("whats the big idea?");

  if (   destination_directory[0] != 0 && 
       ! eventsearch && ! biasout[0] && ! histmode) 
    fprintf(stderr,"%s\n%s\n",
	    "warning: not using the directory if neither",
	    "eventsearch, histmode  nor biasout is selected.");
  if (  destination_directory[0] == 0 && (histmode+eventsearch) > 1 )
    usage("need to specify a directory to access more than one function.");

  if (destination_directory[0]) {
    char *basename;
    /* fix up the output names */
    if (biasout[0]) {
      /* search for directory delimiters */
      if ((basename=strrchr(biasout,(int)'/'))!=NULL) {
	basename++;
	if (strlen(basename)==0) {
	  usage("output name for the bias doesn't make sense.");
	}
      } else basename=biasout;
      sprintf(tempstr,"%s%s",destination_directory,basename);
      sprintf(biasout,"%s",tempstr);
    }
    if (eventsearch || histmode)  {
      /* specify a different filename for each of the eventlists produced */
      for (i=0;i<fni;i++) {
	if ((basename=strrchr(filename[i],(int)'/'))!=NULL) {
	  basename++;
	  if (strlen(basename)==0) {
	    usage("output name for eventfile doesn't make sense.");
	  }
	} else basename=filename[i];

	if (eventsearch) 
	  sprintf(eventfilename[i],
		  "%s%s.evlist",destination_directory,basename);

	if (histmode)
	  sprintf(histfilename[i],
		  "%s%s.pixhist",destination_directory,basename);
      }
    }
  }

  /* generate internal format info. */
  switch (ocspec) {
  case raw:
    break;
  case astd:
    {
      /* set up ocsample -- raw format */
      /*ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=1;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=553-517+1;
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=512; 
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=100;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=516+i+340*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=4+i+340*oc; 
      }
    }
    break;
  case lbox: 
    {
      /* set up ocsample -- currently only LBOX format is supported */
      /*ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=4;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=337-261+1;
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=256; /*can be 256 for faster processing*/
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=100;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=260+i+340*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=4+i+340*oc; 
      }
    }
    break;
  case hirefs: 
    {
      /* set up ocsample -- currently only LBOX format is supported */
      /*ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=2;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=337-261+1;
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=256; /*can be 256 for faster processing*/
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=100;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=260+i+340*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=4+i+340*oc; 
      }
    }
    break;
  case berlin: 
    {
      /* set up ocsample -- currently only LBOX format is supported */
      /*ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=4;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=276-261+1;
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=256; /*can be 256 for faster processing*/
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=100;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=260+i+280*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=4+i+280*oc; 
      }
    }
    break;
  case bnl_e2v_studycontract:
    {
      /* set up ocsample -- 
      /* ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=1;
      //      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=276-261+1;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=2024-2009+1;
      //      ncor[0]=ncor[1]=ncor[2]=ncor[3]=256; /*can be 256 for faster processing*/
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=2008; /*can be 256 for faster processing*/
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=100;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=0+ncor[oc]+i+2048*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=0+i+2048*oc; 
      }
    }
    break;
  case e2v_ccd250:
    {
      /* set up ocsample -- 
      /* ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=1;
      //      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=276-261+1;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=542-(512+10)+1;
      //      ncor[0]=ncor[1]=ncor[2]=ncor[3]=256; /*can be 256 for faster processing*/
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=512+10; /*can be 256 for faster processing*/
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=2002;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=0+ncor[oc]+i+2048*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=0+i+2048*oc; 
      }
    }
    break;
  case lsst_sta_studycontract:
    {
      /* set up ocsample -- 
      /* ocsample[][] contains x-addresses for the OC regions to be sampled.*/
      noc=1;
      //      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=276-261+1;
      nocpix[0]=nocpix[1]=nocpix[2]=nocpix[3]=529-510+1;
      //      ncor[0]=ncor[1]=ncor[2]=ncor[3]=256; /*can be 256 for faster processing*/
      ncor[0]=ncor[1]=ncor[2]=ncor[3]=509; /*can be 256 for faster processing*/
      /* the idea is 100 rows of the overclocks will be used as samples */
      nysample=1999;
      for (i=0;i<nysample;i++)
	ocsample_y[i]=i;
      for (oc=0;oc<noc;oc++) {
	for (i=0;i<nocpix[oc];i++) 
	  ocsample[oc][i]=0+ncor[oc]+i+2048*oc;
	for (i=0;i<ncor[oc];i++)
	  OCcorrect[oc][i]=0+i+2048*oc; 
      }
    }
    break;
  default:
    usage("unspecified frame format.");
  }

  /* first make up lookup arrays for OC sampling and also for correcting */
  for (i=0;i<2048;i++) {    ocsLU[i]=-1;occLU[i]=-1;  }
  for (oc=0;oc<noc;oc++) {
    /* these LUTs allow faster determination of OC parameters. */
    fprintf(stderr,"for oc %d filling in spaces between %d and %d with %d.\n",oc,ocsample[oc][0],ocsample[oc][nocpix[oc]-1],oc);
    fprintf(stderr,"for oc %d filling in spaces between %d and %d with %d.\n",oc,OCcorrect[oc][0],OCcorrect[oc][ncor[oc]-1],oc);
    for (xi=0;xi<nocpix[oc];xi++) ocsLU[ocsample[oc][xi]]=oc;
    for (xi=0;xi<ncor[oc];xi++)   occLU[OCcorrect[oc][xi]]=oc;
    for (fi=0;fi<FNMAX;fi++)      ocval[fi][oc]=0.0;
    for (fi=0;fi<FNMAX;fi++)      oc_int[fi][oc]=0;
  }


  if (biasout[0]) {
    char tmp[2048];
    sprintf(tmp,"%s",biasout);
    sprintf(biasout,"!%s",tmp);
    if (fits_create_file(&ffout,biasout,&status)) 
      printerror(status);
  }
  for (i=1;i<=fni;i++) {
    weight[i]=1.0/(1.0*i);
  }

  /* open up files for output */
  if (eventsearch) {
    for (fi=0;fi<fni;fi++) {
      if (destination_directory[0]) {
	if ((fevlist[fi]=fopen(eventfilename[fi],"w"))==NULL)
	  usage("can't open output evlist file.");
	fclose(fevlist[fi]);
	fevlist[fi]=NULL;
      }	else fevlist[fi]=stdout;
    }
  }

  if (histmode) {
    for (fi=0;fi<fni;fi++) {
      if ((hist[fi]=(int*)malloc((histmax-histmin+1)*sizeof(int)))==NULL)
	usage("can't allocate histogram.");
      for(i=histmin;i<=histmax;i++) hist[fi][i-histmin]=0;
    }
  }
  
  /* now there should be fni files open pointed to by fp[0..fni-1]. */
  //
  //  if (fp[0]==NULL) {
  //    if ((fp[0]=fopen(filename[0],"r"))==NULL) usage("can't open file.");
  //    fseek(fp[0],filepos[0],0);
  //    filepos[0]=0L;
  //  }
  //  read_header(header,fp[0]);
  //  if (filepos[0]=0L) {
  //    filepos[0]=ftell(fp[0]);
  //    fclose(fp[0]);
  //    fp[0]=NULL;
  //  }

  /* immediately output this to the output file. */
  /*  get size info. from the header */
  if (fits_open_file(&ffp[0],filename[0],READONLY,&status)) 
    printerror(status);

  if (biasout[0]) {
    //    if (fits_copy_file(ffp[0],ffout,1,1,1,&status))
    //      printerror(status);
    long naxis[2];
    if (fits_get_img_size(ffp[0],2,naxis,&status)) {
      fprintf(stderr,"can't find NAXIS* in this header unit..\n");
      printerror(status);
    }
    if (fits_create_img(ffout,LONG_IMG,2,naxis,&status)) {
      fprintf(stderr,"can't create image..\n");
      printerror(status);
    }
    // add history string to ffout
    if (fits_write_history(ffout," ",&status))
      printerror(status);
    if (fits_write_history(ffout,
			   "median image created by command:",
			   &status))
      printerror(status);
    if (fits_write_history(ffout,cmd,&status))
      printerror(status);
    if (fits_write_history(ffout," ",&status))
      printerror(status);

  }
  long naxis[2];
  {
    if (fits_get_img_size(ffp[0],2,naxis,&status)) {
      fprintf(stderr,"can't find NAXIS* in this header unit..\n");
      printerror(status);
    }
    nx=naxis[0];    ny=naxis[1];
    fprintf(stderr,"nx = %d ny = %d\n",nx,ny);
  }

  if (input_biasfile[0]) {
    if (fits_get_img_size(ffp[0],2,naxis,&status)) {
      fprintf(stderr,"can't find NAXIS* in this header unit..\n");
      printerror(status);
    }
    if ((naxis[0]-(long)nx != 0) || (naxis[1]-(long)ny != 0)) {
      fprintf(stderr,"naxis[01] from first file: (%d,%d)\n",nx,ny);
      fprintf(stderr,"naxis[01] from input bias: (%ld,%ld)\n",naxis[0],naxis[1]);
      usage("mis-matching header parameters??");
    }
  }

  for(fi=1;fi<fni;fi++)  {

    if (ffp[fi]==NULL) { // if file is not open, open it and set it to close
      if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	printerror(status);
      filepos[fi]=0L;
    }

    if (fits_get_img_size(ffp[fi],2,naxis,&status)) {
      fprintf(stderr,"can't find NAXIS* in this header unit..\n");
      printerror(status);
    }

    if (filepos[fi]==0L) { // close file etc.
      if (fits_close_file(ffp[fni],&status))
	printerror(status);
      ffp[fni]=NULL;
      filepos[fni]=1L;
    }

    if (nx!=naxis[0]) {
      usage("mis-matching header parameters??");
    }
    if (ny!=naxis[1]) {
      usage("mis-matching header parameters??");
    }
  }

  if (burstmode) {
    long fpixel[2];
    fpixel[0]=1L;      fpixel[1]=1L;
    /* in burstmode it should also be possible to find the events */
    memorysize=nx*ny*sizeof(int);
    if ((medianbias=(int*)malloc(memorysize))==NULL)
      usage("can't allocate bias array.");
    if ((occLUTAB=(char*)malloc(nx*ny*sizeof(char)))==NULL)
      usage("can't allocate LUTAB array.");
    /* try something fancy here -- allocate arrays for all of the files */
    for (fi=0;fi<fni;fi++) {
      if ((dp[fi]=(int*)malloc(memorysize))==NULL) 
	usage("can't allocate array.");
    }
    
    /* slurp in the files */
    npix=nx*ny;
    /* first prep the LUTAB */
    for (i=0;i<npix;i++) {
      occLUTAB[i]=(char)(occLU[(i%nx)]+1);
      //      occLUTAB[i]=1; // specify correction for all pixels. (disables nice ocsample stuff)
    }
    fprintf(stderr,"BURSTMODE: reading files..\n");
    if (input_biasfile[0]) {
      if (fits_read_pix(ffib,TINT,fpixel,(long)nx*ny,NULL,medianbias,
			NULL,&status)) {
	fprintf(stderr,"can't load the specified bias.");
	printerror(status);
      }
    }
    fprintf(stderr,"will read in %d files..\n",fni);
    for (fi=0;fi<fni;fi++) {

      if (ffp[fi]==NULL) {
	if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	  printerror(status);
	fprintf(stderr,"just opened file %s.. ",filename[fi]);
	filepos[fi]=0L;
      }
      if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx*ny,NULL,dp[fi],
			NULL,&status)) {
	fprintf(stderr,"can't load the specified file..");
	printerror(status);
      }
      if (filepos[fi]==0L) {
	if (fits_close_file(ffp[fni],&status))
	  printerror(status);
	fprintf(stderr,"and closed it.\n",filename[fi]);
	ffp[fni]=NULL;
	filepos[fni]=1L;
      }
      //      fprintf(stderr,"done.\n");
    }

    if (OCcorrection) {
      /* before proceeding, look into the overclocks and apply 
	 corrections in place. */
      fprintf(stderr,"doing something with OC correction..\n");
      evaluate_OC_vals(dp,fni,nx,ny,nysample,ocsample_y,noc,
		       nocpix,ocsample,occLUTAB,ocval);
      for (fi=0;fi<fni;fi++)
	for (oc=0;oc<noc;oc++)
	  oc_int[fi][oc]=(int)(floor(ocval[fi][oc]+0.5));
    } 
    
    if (input_biasfile[0]) {
      for (i=0;i<npix;i++) {
	if (occLUTAB[i]) {
	  for(fi=0;fi<fni;fi++) {
	    dp[fi][i] -= oc_int[fi][occLUTAB[i]-1];
	    pix[fi]=dp[fi][i];
	  }
	  median(pix,fni,&tmp);
	  for(fi=0;fi<fni;fi++)
	    dp[fi][i]-=medianbias[i];
	  medianbias[i]=tmp-medianbias[i];
	}
      }
    } else {
      /* calculate the bias frame, and output the difference.
	 use the input_bias for subsequent analysis. */
      for (i=0;i<npix;i++) {
	tmp=0;
	if (occLUTAB[i]) {
	  for(fi=0;fi<fni;fi++) {
	    dp[fi][i]-=oc_int[fi][occLUTAB[i]-1];
	    pix[fi]=dp[fi][i];
	  }
	  median(pix,fni,&tmp);
	  for(fi=0;fi<fni;fi++)
	    dp[fi][i]-=tmp;
	}
	medianbias[i]=tmp; 
      }
      
      if (biasout[0]) {
	long fpixel[2],lpixel[2];
	
	fpixel[0]=1L;	        fpixel[1]=1L;
	lpixel[0]=(long)nx;	lpixel[1]=(long)ny;
	
	if (fits_write_subset(ffout,TINT,fpixel,lpixel,medianbias,&status))
	  printerror(status);
      }
    }
    if (0) {
    if (input_biasfile[0]) {
      for (i=0;i<npix;i++) {
	if (occLUTAB[i]) {
	  for(fi=0;fi<fni;fi++) 
	    dp[fi][i]-=(medianbias[i]+oc_int[fi][occLUTAB[i]-1]);
	}
      }
    } else {
      for (i=0;i<npix;i++) {
	if (occLUTAB[i]) {
	  for(fi=0;fi<fni;fi++) {
	    dp[fi][i]-=oc_int[fi][occLUTAB[i]-1];
	    pix[fi]=dp[fi][i];
	  }
	  median(pix,fni,&tmp);
	  for(fi=0;fi<fni;fi++)
	    dp[fi][i]-=tmp;
	}
      }
    }
    }

    if (medianbias) free(medianbias);

    /* now it should be possible to just analyse the *dp pointers */
    if (histmode) {
      for(fi=0;fi<fni;fi++) {
	if (destination_directory[0]) {
	  if ((fhist[fi]=fopen(histfilename[fi],"w"))==NULL)
	    usage("can't open output histogram file.");
	} else {
	  fhist[fi]=stdout;
	}
	for (i=0;i<npix;i++) {
	  if (occLUTAB[i]) {
	    if (((tmp=dp[fi][i])>=histmin) && (tmp<=histmax)) {
	      hist[fi][tmp-histmin]++;
	    }
	  }
	}
	for (i=histmin;i<=histmax;i++) {
	  if (i==0) continue;
	  fprintf(fhist[fi],"%d %d\n",i,hist[fi][i-histmin]);
	}
	if (destination_directory[0]) {
	  if (fhist[fi] != stdout) {
	    fclose(fhist[fi]);
	  }
	}
      }
    }

    
    if (eventsearch) {
      for(fi=0;fi<fni;fi++) {
	//	fprintf(stderr,"going through array for file %s..\n",filename[fi]);
	nev=0;
	/* examine contents of the median-subtracted frame for events. */
	/* examine 3 rows of the image file at a time. */
	if (1) { // rebin in place by 2
	  int *rebin_data=NULL;
	  int rnx,rny;
	  rnx=ceil(nx/(float)REB);
	  rny=ceil(ny/(float)REB);
	  if ((rebin_data=(int*)malloc(rnx*rny*sizeof(int)))==NULL) {
	    usage("can't allocate rebin_data!\n");
	  }
	  i=rnx*rny;	  while (i--) rebin_data[i]=0;
	  {
	    int k,j;
	    for (k=0;k<nx;k++) {
	      for (j=0;j<ny;j++) {
		if (occLUTAB[k+j*nx]) 
		  rebin_data[(int)floor(k/(float)REB)+(int)floor(j/(float)REB)*rnx] += dp[fi][k+j*nx];
	      }
	    }
	    for (j=1;j<rny-1;j++) {
	      for (k=1;k<rnx-1;k++) {
		i=k+j*rnx;
		cpix=&(rebin_data[i]);
		if (*cpix < (int) evthresh) continue;
		if (*cpix >= *(cpix+rnx)   &&     *cpix >= *(cpix+1)    && 
		    *cpix >  *(cpix-1)     &&     *cpix >  *(cpix-rnx)   &&
		    *cpix >= *(cpix+rnx+1) &&     *cpix >= *(cpix+rnx-1) &&
		    *cpix >  *(cpix-rnx+1) &&     *cpix >  *(cpix-rnx-1)) {
		  /* local maximum hit -- output an evlist ! */
		  x=i%rnx;y=i/rnx;
		  nev++; event.x=x;event.y=y;
		  for (yi=-1;yi<=1;yi++)
		    for (xi=-1;xi<=1;xi++) 
		      event.data[xi+1+(yi+1)*3]=(int)(*(cpix+yi*rnx+xi));
		  pushevent(&event,&(evstack[fi]));
		}
	      }
	    }
	  }
	  free(rebin_data);
	} else {
	  for (i=nx+1;i<npix-nx-1;i++) {
	    if (!occLUTAB[i]) continue;
	    cpix = (int*)dp[fi] + i;
	    if (*cpix < (int) evthresh) continue;
	    if (*cpix >= *(cpix+nx)   &&     *cpix >= *(cpix+1)    && 
		*cpix >  *(cpix-1)    &&     *cpix >  *(cpix-nx)   &&
		*cpix >= *(cpix+nx+1) &&     *cpix >= *(cpix+nx-1) &&
		*cpix >  *(cpix-nx+1) &&     *cpix >  *(cpix-nx-1)) {
	      /* local maximum hit -- output an evlist ! */
	      x=i%nx;y=i/nx;
	      nev++; event.x=x;event.y=y;
	      for (yi=-1;yi<=1;yi++)
		for (xi=-1;xi<=1;xi++) 
		  event.data[xi+1+(yi+1)*3]=(int)(*(cpix+yi*nx+xi));
	      pushevent(&event,&(evstack[fi]));
	      i++; 
	    }
	  }
	}
	fprintf(stderr," got %d events.\n",nev);
	/* and close the eventfile if thats where output has been directed. */
      }
    }
  } else {     
    /* should do this 3 lines at a time. that way, eventsearch
       and histmode will be possible. probably faster all around too. 
    */
    int *rp[FNMAX][3],*median_row[3],*tpix,*cpix,*bpix,tmp_med_row[2048];
    int topindex,cenindex,botindex,row;
    long fpixel[2];

    fprintf(stderr,"this is not burst mode!!\n");
    if (input_biasfile[0] && !biasout[0]) {
      /* can do this faster if there will never be a need to calculate median. 
	 also, can do one file at a time -- no limit to number of files 
	 processed. */

      for (fi=0;fi<fni;fi++) 
	/* allocate pointers to rows. 3 rows per file */
	for (row=0;row<3;row++) {
	  if ((rp[fi][row]=(int*)malloc(nx*sizeof(int)))==NULL)
	    usage("can't allocate row array??");
	}
	   
      if (OCcorrection) {
	/* first find the overclock values */
	evaluate_file_OC_vals(ffp,fni,nx,ny,nysample,ocsample_y,
			      noc,nocpix,ocsample,occLUTAB,ocval);
	for (fi=0;fi<fni;fi++)
	  for (oc=0;oc<noc;oc++)
	    oc_int[fi][oc]=(int)floor(ocval[fi][oc]);
      }
      
      npix=0;
	   
      /* first read in 2 rows for each file. and median the pixels */
      for (row=0;row<2;row++) {
	fpixel[0]=1L;      fpixel[1]=(long)(row+1);
	for (fi=0;fi<fni;fi++) {

	  if (ffp[fi]==NULL) {
	    if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	      printerror(status);
	    filepos[fi]=0L;
	  }
	  if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx,NULL,rp[fi][row],
			    NULL,&status)) {
	    fprintf(stderr,"can't load the specified file..");
	    printerror(status);
	  }
	  if (filepos[fi]==0L) {
	    if (fits_close_file(ffp[fni],&status))
	      printerror(status);
	    ffp[fni]=NULL;
	    filepos[fni]=1L;
	  }
	}
	     
	if (fits_read_pix(ffib,TINT,fpixel,(long)nx,NULL,tmp_med_row,
			  NULL,&status)) {
	  fprintf(stderr,"can't load the specified bias.");
	  printerror(status);
	}
	
	for (i=0;i<nx;i++) {
	  if (((oc=occLU[i])+1)) {
	    for (fi=0;fi<fni;fi++) {
	      rp[fi][row][i]-=oc_int[fi][oc];
	      rp[fi][row][i]-=tmp_med_row[i]; /* subtract input */
	    }
	  } else 
	    for (fi=0;fi<fni;fi++)
	      rp[fi][row][i]=0;
	}
	
	if (histmode) {
	  for (i=0;i<nx;i++)
	    for (fi=0;fi<fni;fi++)
	      if ((occLU[i]+1)&&((tmp=rp[fi][row][i])>=histmin)&&(tmp<=histmax)) 
		hist[fi][tmp-histmin]++;
	}
      }
      
      /* now there exist medians for the first two data rows. 
	 from here on need to perform a caterpillar type analysis through
	 the frames. */
	   
      for (row=2;row<ny;row++) {
	     
	fpixel[0]=1L;      fpixel[1]=(long)(row+1);
	topindex=row%3;
	cenindex=(row-1)%3;
	botindex=(row-2)%3;
	
	for (fi=0;fi<fni;fi++) {

	  if (ffp[fi]==NULL) {
	    if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	      printerror(status);
	    filepos[fi]=0L;
	  }
	  if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx,NULL,rp[fi][topindex],
			    NULL,&status)) {
	    fprintf(stderr,"can't load the specified file..");
	    printerror(status);
	  }
	  if (filepos[fi]==0L) {
	    if (fits_close_file(ffp[fni],&status))
	      printerror(status);
	    ffp[fni]=NULL;
	    filepos[fni]=1L;
	  }
	  
	}
	
	//	  usage("can't read median file??");
	//
	// read a row of the input bias file
	//

	if (fits_read_pix(ffib,TINT,fpixel,(long)nx,NULL,tmp_med_row,
			  NULL,&status)) {
	  fprintf(stderr,"can't load the specified bias.");
	  printerror(status);
	}
	
	for (i=0;i<nx;i++) {
	  if (((oc=occLU[i])+1)) {
	    for (fi=0;fi<fni;fi++) {
	      rp[fi][topindex][i]-=oc_int[fi][oc];
	      rp[fi][topindex][i]-=tmp_med_row[i]; /* subtract input */
	    }
	  } else 
	    for (fi=0;fi<fni;fi++)
	      rp[fi][topindex][i]=0;
	}
	
	if (histmode) {
	  for (i=0;i<nx;i++)
	    for (fi=0;fi<fni;fi++)
	      if ((occLU[i]+1)&&
		  ((tmp=rp[fi][topindex][i])>=histmin)&&
		  (tmp<=histmax)) 
		hist[fi][tmp-histmin]++;
	}
	
	if (! (eventsearch || histmode)) continue;  /* dont waste time */
	
	/* now have 3 rows to deal with. 
	   median biasses have already been subtracted.
	   top    row is rp[fi][topindex]
	   center row is rp[fi][cenindex]
	   bottom row is rp[fi][botindex]
	*/
	
	if (eventsearch) {
	  for (fi=0;fi<fni;fi++) {
	    tpix=&(rp[fi][topindex][0]);
	    cpix=&(rp[fi][cenindex][0]);
	    bpix=&(rp[fi][botindex][0]);
	    for (x=1;x<nx-1;x++) {
	      tpix++;cpix++;bpix++;
	      if (*(cpix)>=evthresh  &&
		  *(cpix)>=*(cpix+1) && *(cpix)>=*(tpix)   &&
		  *(cpix)> *(cpix-1) && *(cpix)> *(bpix)   &&
		  *(cpix)>=*(tpix+1) && *(cpix)>=*(tpix-1) &&
		  *(cpix)> *(bpix+1) && *(cpix)> *(bpix-1)) {
		/* output the event. */
		event.x=x;event.y=row;
		for (xi=-1;xi<=1;xi++) {
		  event.data[1+xi]=(int)(*(bpix+xi));
		  event.data[4+xi]=(int)(*(cpix+xi));
		  event.data[7+xi]=(int)(*(tpix+xi));
		}
		event.framenum=fi;
		event.chipnum=0;
		pushevent(&event,&(evstack[fi]));
	      }
	    }
	  }
	}
      }
	   
      /* flush all pending histogram data here. */
      if (histmode) {
	for (fi=0;fi<fni;fi++) {
	  if (destination_directory[0]) {
	    if ((fhist[fi]=fopen(histfilename[fi],"w"))==NULL)
	      usage("can't open output histogram file.");
	  }	else fhist[fi]=stdout;
	  
	  for (i=histmin;i<=histmax;i++) {
	    if (i==0) continue;
	    fprintf(fhist[fi],"%d %d\n",i,hist[fi][i-histmin]);
	  }
	  if (destination_directory[0]) {
	    fclose(fhist[fi]);
	  }
	}
      }
      
      /* release the space */
      for (row=0;row<3;row++) {
	for (fi=0;fi<fni;fi++)
	  if (rp[fi][row]) 
	    free(rp[fi][row]);
      }
      
    } else { /* ! (input_biasfile[0] && !biasout[0])  */

      fprintf(stderr,"so far so good (1)\n");

      for (fi=0;fi<fni;fi++) 
	/* allocate pointers to rows. 3 rows per file */
	for (row=0;row<3;row++) {
	  if ((rp[fi][row]=(int*)malloc(nx*sizeof(int)))==NULL)
	    usage("can't allocate row array??");
	  if ((median_row[row]=(int*)malloc(nx*sizeof(int)))==NULL)
	    usage("can't allocate median row array??");
	  /* initialize median_row */
	  for (i=0;i<nx;i++)
	    median_row[row][i]=0;
	}
      
      fprintf(stderr,"so far so good (2)\n");

      if (OCcorrection) {
	/* first find the overclock values */
	evaluate_file_OC_vals(ffp,fni,nx,ny,nysample,ocsample_y,
			      noc,nocpix,ocsample,occLUTAB,ocval);
	for (fi=0;fi<fni;fi++)
	  for (oc=0;oc<noc;oc++)
	    oc_int[fi][oc]=(int)floor(ocval[fi][oc]);
      }
      
      fprintf(stderr,"so far so good (3)\n");

      npix=0;
      
      /* first read in 2 rows for each file. and median the pixels */
      for (row=0;row<2;row++) {
	fpixel[0]=1L;      fpixel[1]=(long)(row+1);
	for (fi=0;fi<fni;fi++) {
	  
	  //	  if (fp[fi]==NULL) {
	  //	    if ((fp[fi]=fopen(filename[fi],"r"))==NULL) 
	  //	      usage("can't open file.");
	  //	    fseek(fp[fi],filepos[fi],0);
	  //	    filepos[fi]=0L;
	  //	  }
	  //	  
	  //	  if (fread(rp[fi][row],sizeof(int),nx,fp[fi])!=nx)
	  //	    usage("can't read the first two rows??");
	  //	  
	  //	  if (filepos[fi]==0L) {
	  //	    filepos[fi]=ftell(fp[fi]);
	  //	    fclose(fp[fi]);
	  //	    fp[fi]=NULL;
	  //	  }
	  //	  
	
	  fprintf(stderr,"so far so good (row=%d,fi=%d)\n",row,fi);

	  if (ffp[fi]==NULL) {
	    if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	      printerror(status);
	    filepos[fi]=0L;
	  }
	  if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx,NULL,rp[fi][row],
			    NULL,&status)) {
	    fprintf(stderr,"can't load the specified file..");
	    printerror(status);
	  }
	  if (filepos[fi]==0L) {
	    if (fits_close_file(ffp[fni],&status))
	      printerror(status);
	    ffp[fni]=NULL;
	    filepos[fni]=1L;
	  }
	}

	fprintf(stderr,"so far so good (4)\n",row,fi);

	if (input_biasfile[0]) {
	  //	  if (fread(tmp_med_row,sizeof(int),nx,fib)!=nx)
	  //	    usage("can't read median file??");
	  if (fits_read_pix(ffib,TINT,fpixel,(long)nx,NULL,tmp_med_row,
			    NULL,&status)) {
	    fprintf(stderr,"can't load the specified bias.");
	    printerror(status);
	  }
	  
	  for (i=0;i<nx;i++) {
	    if (((oc=occLU[i])+1)) {
	      for (fi=0;fi<fni;fi++) {
		rp[fi][row][i]-=oc_int[fi][oc];
		
		pix[fi]=rp[fi][row][i];
	      }
	      if (biasout[0]) {
		median(pix,fni,median_row[row]+i);
		median_row[row][i]-=tmp_med_row[i];
	      }
	      for (fi=0;fi<fni;fi++)
		rp[fi][row][i]-=tmp_med_row[i]; /* subtract input */
	    } else 
	      for (fi=0;fi<fni;fi++)
		rp[fi][row][i]=0;
	  }
	  
	} else {
	  fprintf(stderr,"so far so good (5)\n");
	  /* usual on-the-fly bias determination */
	       for (i=0;i<nx;i++) {
		 if (((oc=occLU[i])+1)) {
		   for (fi=0;fi<fni;fi++) {
		     rp[fi][row][i]-=oc_int[fi][oc];
		     pix[fi]=rp[fi][row][i];
		   }
		   
		   fprintf(stderr,"so far so good (will find median of row %d\n",row);
		   median(pix,fni,median_row[row]+i);
		   //		   fprintf(stderr,"median pixel turns out: %d\n",*(median_row[row]+i));
		   for (fi=0;fi<fni;fi++) {
		     rp[fi][row][i]-=median_row[row][i];
		   }
		 } else 
		   //		   fprintf(stderr,"didn't compute median for this..\n");
		   for (fi=0;fi<fni;fi++)
		     rp[fi][row][i]=0;
	       }
	     }
	     
	if (biasout[0]) {
	  long fpixel[2],lpixel[2];

	  fpixel[0]=1L;	        fpixel[1]=(long)(row+1);
	  lpixel[0]=(long)nx;	lpixel[1]=(long)(row+1);

	  if (fits_write_subset(ffout,TINT,
				fpixel,lpixel,median_row[row],&status))
	    printerror(status);
	}
	  

	     if (histmode) {
	       for (i=0;i<nx;i++)
		 for (fi=0;fi<fni;fi++)
		   if ((occLU[i]+1)&&((tmp=rp[fi][row][i])>=histmin)&&(tmp<=histmax)) 
		     hist[fi][tmp-histmin]++;
	     }
	   }
	   
	   /* now there exist medians for the first two data rows. 
	      from here on need to perform a caterpillar type analysis through
	      the frames. */
	   
      //      fprintf(stderr,"not burst mode.. 7\n");

      for (row=2;row<ny;row++) {
	fpixel[0]=1L;      fpixel[1]=(long)(row+1);
	
	topindex=row%3;
	cenindex=(row-1)%3;
	botindex=(row-2)%3;
	
	for (fi=0;fi<fni;fi++) {
	  
	  //	       if (fp[fi]==NULL) {
	  //		 if ((fp[fi]=fopen(filename[fi],"r"))==NULL) 
	  //		   usage("can't open file.");
	  //		 fseek(fp[fi],filepos[fi],0);
	  //		 filepos[fi]=0L;
	  //	       }
	  //
	  //	       if (fread(rp[fi][topindex],sizeof(int),nx,fp[fi])!=nx)
	  //		 usage("can't read in specified row??");
	  //
	  //	       if (filepos[fi]==0L) {
	  //		 filepos[fi]=ftell(fp[fi]);
	  //		 fclose(fp[fi]);
	  //		 fp[fi]=NULL;
	  //	       }
	  if (ffp[fi]==NULL) {
	    if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	      printerror(status);
	    filepos[fi]=0L;
	  }
	  if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx,NULL,rp[fi][topindex],
			    NULL,&status)) {
	    fprintf(stderr,"can't load the specified file..");
	    printerror(status);
	  }
	  if (filepos[fi]==0L) {
	    if (fits_close_file(ffp[fni],&status))
	      printerror(status);
	    ffp[fni]=NULL;
	    filepos[fni]=1L;
	  }
	}
	
	if (input_biasfile[0]) {
	  //	       if (fread(tmp_med_row,sizeof(int),nx,fib)!=nx)
	  //		 usage("can't read median file??");
	  if (fits_read_pix(ffib,TINT,fpixel,(long)nx,NULL,tmp_med_row,
			    NULL,&status)) {
	    fprintf(stderr,"can't load the specified bias.");
	    printerror(status);
	  }
	  
	  for (i=0;i<nx;i++) {
	    if (((oc=occLU[i])+1)) {
	      for (fi=0;fi<fni;fi++) {
		rp[fi][topindex][i]-=oc_int[fi][oc];
		pix[fi]=rp[fi][topindex][i];
	      }
	      if (biasout[0]) {
		median(pix,fni,median_row[topindex]+i);
		median_row[topindex][i]-=tmp_med_row[i];
	      }
	      for (fi=0;fi<fni;fi++) 
		rp[fi][topindex][i]-=tmp_med_row[i]; /* subtract input */
	    } else 
	      for (fi=0;fi<fni;fi++)
		rp[fi][topindex][i]=0;
	  }
	  
	} else {
	  /* usual on-the-fly bias determination */
	  for (i=0;i<nx;i++) {
	    if (((oc=occLU[i])+1)) {
	      for (fi=0;fi<fni;fi++) {
		rp[fi][topindex][i]-=oc_int[fi][oc];
		pix[fi]=rp[fi][topindex][i];
	      }
	      
	      median(pix,fni,median_row[topindex]+i);
	      for (fi=0;fi<fni;fi++) {
		rp[fi][topindex][i]-=median_row[topindex][i];
	      }
	    } else 
	      for (fi=0;fi<fni;fi++)
		rp[fi][topindex][i]=0;
	  }
	}
	
	if (biasout[0]) { // output the row
	  long fpixel[2],lpixel[2];
	  
	  fpixel[0]=1L;	        fpixel[1]=(long)(row+1);
	  lpixel[0]=(long)nx;	lpixel[1]=(long)(row+1);
	  
	  if (fits_write_subset(ffout,TINT,
				fpixel,lpixel,median_row[topindex],&status))
	    printerror(status);
	}
	
	if (histmode) {
	  for (i=0;i<nx;i++)
	    for (fi=0;fi<fni;fi++)
	      if ((occLU[i]+1)&&
		  ((tmp=rp[fi][topindex][i])>=histmin)&&
		  (tmp<=histmax)) 
		hist[fi][tmp-histmin]++;
	}
	
	if (! (eventsearch || histmode)) continue;  /* dont waste time */
	
	/* now have 3 rows to deal with. 
	   median biasses have already been subtracted.
	   top    row is rp[fi][topindex]
	   center row is rp[fi][cenindex]
	   bottom row is rp[fi][botindex]
	*/
	
	if (eventsearch) {
	  for (fi=0;fi<fni;fi++) {
	    tpix=&(rp[fi][topindex][0]);
	    cpix=&(rp[fi][cenindex][0]);
	    bpix=&(rp[fi][botindex][0]);
	    for (x=1;x<nx-1;x++) {
	      tpix++;cpix++;bpix++;
	      if (*(cpix)>=evthresh  &&
		  *(cpix)>=*(cpix+1) && *(cpix)>=*(tpix)   &&
		  *(cpix)> *(cpix-1) && *(cpix)> *(bpix)   &&
		  *(cpix)>=*(tpix+1) && *(cpix)>=*(tpix-1) &&
		  *(cpix)> *(bpix+1) && *(cpix)> *(bpix-1)) {
		/* output the event. */
		event.x=x;event.y=row;
		for (xi=-1;xi<=1;xi++) {
		  event.data[1+xi]=(int)(*(bpix+xi));
		  event.data[4+xi]=(int)(*(cpix+xi));
		  event.data[7+xi]=(int)(*(tpix+xi));
		}
		event.framenum=fi;
		event.chipnum=0;
		//		fprintf(stderr,"pushing event.. fi = %d\n",fi);
		pushevent(&event,&(evstack[fi]));
	      }
	    }
	  }
	}
      }
      //      fprintf(stderr,"not burst mode.. 8\n");

	   
      /* flush all pending histogram data here. */
      if (histmode) {
	for (fi=0;fi<fni;fi++) {
	  if (destination_directory[0]) {
	    if ((fhist[fi]=fopen(histfilename[fi],"w"))==NULL)
	      usage("can't open output histogram file.");
	  }	else fhist[fi]=stdout;
	  
	  for (i=histmin;i<=histmax;i++)
	    fprintf(fhist[fi],"%d %d\n",i,hist[fi][i-histmin]);
	  
	  if (destination_directory[0]) {
	    fclose(fhist[fi]);
	  }
	}
      }

      //      fprintf(stderr,"not burst mode.. 9\n");
      
      /* release the space */
      for (row=0;row<3;row++) {
	if (median_row[row])
	  free(median_row[row]);
	for (fi=0;fi<fni;fi++)
	  if (rp[fi][row]) 
	    free(rp[fi][row]);
      }
      //      fprintf(stderr,"not burst mode.. 10\n");
    }  
  } /* end of normal (non-burst) mode */

  //  fprintf(stderr,"?? burst mode.. 11\n");
  /* spray out evlists into appropriate channels */
  if (eventsearch) {
    if (destination_directory[0]) {
      //      fprintf(stderr,"?? burst mode.. 11a\n");
      for (fi=0;fi<fni;fi++) {
	if ((fevlist[fi]=fopen(eventfilename[fi],"w"))==NULL)
	  usage("can't open output evlist file.");
	/* and do the stuff with the evstack */
	current_evstack_pointer=evstack[fi];
	while (current_evstack_pointer!=NULL) {
	  fwrite(
		 &(current_evstack_pointer->event),
		 datastr_size,1,fevlist[fi]);
	  current_evstack_pointer=current_evstack_pointer->prev_ev;
	}
	pulldown_events(evstack[fi]);
	fclose(fevlist[fi]);
      }
    } else {
      //      fprintf(stderr,"?? burst mode.. 11b\n");
      for (fi=0;fi<fni;fi++) {
	// send these all to the stdout..
	/* and do the stuff with the evstack */
	current_evstack_pointer=evstack[fi];
	while (current_evstack_pointer!=NULL) {
	  fwrite(
		 &(current_evstack_pointer->event),
		 datastr_size,1,fevlist[fi]);
	  current_evstack_pointer=current_evstack_pointer->prev_ev;
	}
	pulldown_events(evstack[fi]);
      }
    }
  }
  //  fprintf(stderr,"?? burst mode.. 12\n");

  if (biasout[0]) 
    if (fits_close_file(ffout,&status))
      printerror(status);

  if (burstmode)
    for (fi=0;fi<fni;fi++) if (dp[fi]) free(dp[fi]);
  if (histmode)
    for (fi=0;fi<fni;fi++) if (hist[fi]) free(hist[fi]);
  if (occLUTAB) 
    free (occLUTAB);
}

void median(int x[],int n,int *xmed) {
	int n2,n2p;
	void sort_uint();

	if (n==1) {
	  *xmed=x[0];
	  return;
	}
	//	sort_uint(n,x);
	qsort((void*)x,(size_t)n,sizeof(int),comp_int);
	n2p=(n2=n/2)+1;
	*xmed=(n % 2 ? x[n2p] : (int)floor(0.5*(x[n2]+x[n2p])));
}

int comp_int(const void *a,const void *b) {
  if (*((int*)a) == *((int*)b)) return(0);
  if (*((int*)a) > *((int*)b)) {
    return(1);
  } else {
    return(-1);
  }
}

evaluate_OC_vals 
  (int *dp[],int fni,int nx,int ny,int nysample,int ocsample_y[],int noc,
   int nocpix[],int ocsample[][1024],char occLUTAB[],float ocval[][OCMAX]   
   /*other than the modified *dp arrays, *ocval[] is the return val.*/
   )
//int *dp[];
//int fni,nx,ny;
//int nysample,ocsample_y[],noc,nocpix[],ocsample[][1024];
//char occLUTAB[];
//float ocval[][OCMAX];
{
  int fi,oc,i,yi,xi;
  int min[OCMAX],ochists[OCMAX][OCHISTMAX],psum[OCMAX],nsum[OCMAX],tx,ty,val,
  npix,index;
  int tmp,ocint[FNMAX][OCMAX];

  fprintf(stderr,"sampling mean OC values:\n");
  for (fi=0;fi<fni;fi++) {
    
    /* get an idea for the distribution centers using the 50th row */

    ty=50;
    for (oc=0;oc<noc;oc++) { psum[oc]=0; nsum[oc]=0;  }
    for (oc=0;oc<noc;oc++) {
      index=ty*nx+ocsample[oc][0];
      for (xi=0;xi<nocpix[oc];xi++) {
	psum[oc]+=(int)(dp[fi][index]);
	nsum[oc]++;
	index++;
      }
    }
    
    for (oc=0;oc<noc;oc++) 
      min[oc]=(int)(floor(psum[oc]/(1.0*nsum[oc]))-OCHISTOS); 
    
    /* now back up and make entire histograms for the OC regions */
    /* initialize.. */
    for (oc=0;oc<noc;oc++)
      for (i=0;i<OCHISTMAX;i++) 
	ochists[oc][i]=0;
    /* and sample */
    for (yi=0;yi<nysample;yi++) {
      ty=ocsample_y[yi];
      for (oc=0;oc<noc;oc++) {
	index=ty*nx+ocsample[oc][0];
	for (xi=0;xi<nocpix[oc];xi++) {
	  val=(int)(dp[fi][index])-min[oc];
	  if (val>=0 && val<OCHISTMAX)
	    ochists[oc][val]++;
	  index++;
	}
      }
    }
    
    /* now evaluate the histograms to get OC corrections */
    for (oc=0;oc<noc;oc++) {
      psum[oc]=0;nsum[oc]=0;
      for (i=0;i<OCHISTMAX;i++) {
	psum[oc]+=(i+min[oc])*ochists[oc][i];
	nsum[oc]+=ochists[oc][i];
      }
      ocval[fi][oc]=psum[oc]/(1.0*nsum[oc]);
      ocint[fi][oc]=(int)floor(ocval[fi][oc]+0.5);
    }
    
    /* and here correct for OC fluctuations. */
  }
  
  /* determine variances in OC levels */
  for (oc=0;oc<noc;oc++) {
    int n;
    float mean,variance;
    n=0;
    mean=0.0;
    variance=0.0;
    for (fi=0;fi<fni;fi++) {
      mean+=ocval[fi][oc];
      n++;
    }
    mean /= (float)(n);
    for (fi=0;fi<fni;fi++) {
      fprintf(stderr,"file %d: ocval[%d] = %f (%s)\n",fi,oc,ocval[fi][oc],filename[fi]);
      variance+=pow(ocval[fi][oc]-mean,2.0);
      n++;
    }
    variance = sqrt(variance/(float)(n));
    fprintf(stderr,"OC %d: mean %f variance %f\n",oc,mean,variance);
  }
}


evaluate_file_OC_vals
  (fitsfile *ffp[],int fni,int nx,int ny,int nysample,int ocsample_y[],int noc,
   int nocpix[],int ocsample[][1024],char occLUTAB[],float ocval[][OCMAX]   
   /*other than the modified *dp arrays, *ocval[] is the return val.*/
   )
//fitsfile *ffp[];
//int fni,nx,ny;
//int nysample,ocsample_y[],noc,nocpix[],ocsample[][1024];
//char occLUTAB[];
//float ocval[][OCMAX];
{
  int fi,oc,i,yi,xi;
  int status=0;
  int min[OCMAX],ochists[OCMAX][OCHISTMAX],psum[OCMAX],nsum[OCMAX],tx,ty,val,
  npix,index;
  int tmp,ocint[FNMAX][OCMAX],line[2048];
  long  savepos[FNMAX];
  char str[1024];
  long fpixel[2];

  fprintf(stderr,"sampling mean OC values:\n");

  for (fi=0;fi<fni;fi++) {
    
    //    if (fp[fi]==NULL) {
    //      if ((fp[fi]=fopen(filename[fi],"r"))==NULL) 
    //	usage("can't open file.");
    //      fseek(fp[fi],filepos[fi],0);
    //      filepos[fi]=0L;
    //    }
    
    //    savepos[fi]=ftell(fp[fi]);
    //
    //    fseek(fp[fi],(long)(savepos[fi]+ty*nx*sizeof(int)),0);
    //    fread(line,sizeof(int),nx,fp[fi]);

    //    /* get an idea for the distribution centers using the 50th row */

    ty=50;

    fpixel[0]=1L;
    fpixel[1]=(long)(ty+1);

    if (ffp[fi]==NULL) {
      if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	printerror(status);
      filepos[fi]=0L;
    }
    if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx,NULL,line,
		      NULL,&status)) {
      fprintf(stderr,"can't load the specified file..");
      printerror(status);
    }
    if (filepos[fi]==0L) {
      if (fits_close_file(ffp[fni],&status))
	printerror(status);
      ffp[fni]=NULL;
      filepos[fni]=1L;
    }
    
    for (oc=0;oc<noc;oc++) { psum[oc]=0; nsum[oc]=0;  }
    for (oc=0;oc<noc;oc++) {
      index=ocsample[oc][0];
      for (xi=0;xi<nocpix[oc];xi++) {
	psum[oc]+=(int)line[index];
	nsum[oc]++;
	index++;
      }
      min[oc]=(int)(floor(psum[oc]/(1.0*nsum[oc]))-OCHISTOS); 
    }

    /* now back up and make entire histograms for the OC regions */
    /* initialize.. */
    //    fseek(fp[fi],savepos[fi],0);

    for (oc=0;oc<noc;oc++)
      for (i=0;i<OCHISTMAX;i++) 
	ochists[oc][i]=0;
    /* and sample */
    for (yi=0;yi<nysample;yi++) {
      ty=ocsample_y[yi];

      fpixel[0]=1L;
      fpixel[1]=(long)(ty+1);
      //      fread(line,sizeof(int),nx,fp[fi]);
      if (ffp[fi]==NULL) {
	if (fits_open_file(&ffp[fi],filename[fi],READONLY,&status))
	  printerror(status);
	filepos[fi]=0L;
      }
      if (fits_read_pix(ffp[fi],TINT,fpixel,(long)nx,NULL,line,
			NULL,&status)) {
	fprintf(stderr,"can't load the specified file..");
	printerror(status);
      }
      if (filepos[fi]==0L) {
	if (fits_close_file(ffp[fni],&status))
	  printerror(status);
	ffp[fni]=NULL;
	filepos[fni]=1L;
      }

      for (oc=0;oc<noc;oc++) {
	index=ocsample[oc][0];
	for (xi=0;xi<nocpix[oc];xi++) {
	  val=(int)line[index]-min[oc];
	  if (val>=0 && val<OCHISTMAX)
	    ochists[oc][val]++;
	  index++;
	}
      }
    }

    /* now evaluate the histograms to get OC corrections */
    for (oc=0;oc<noc;oc++) {
      psum[oc]=0;nsum[oc]=0;
      for (i=0;i<OCHISTMAX;i++) {
	psum[oc]+=(i+min[oc])*ochists[oc][i];
	nsum[oc]+=ochists[oc][i];
      }
      ocval[fi][oc]=psum[oc]/(1.0*nsum[oc]);
      ocint[fi][oc]=(int)floor(ocval[fi][oc]+0.5);
    }
    
    /* and here correct for OC fluctuations. */
    //    fseek(fp[fi],savepos[fi],0);
    //    
    //    if (filepos[fi]==0L) {
    //      filepos[fi]=ftell(fp[fi]);
    //      fclose(fp[fi]);
    //      fp[fi]=NULL;
    //    }
    
  }
  
  /* determine variances in OC levels */
  for (oc=0;oc<noc;oc++) {
    int n;
    float mean,variance;
    n=0;
    mean=0.0;
    variance=0.0;
    for (fi=0;fi<fni;fi++) {
      mean+=ocval[fi][oc];
      n++;
    }
    mean /= (float)(n);
    for (fi=0;fi<fni;fi++) {
      variance+=pow(ocval[fi][oc]-mean,2.0);
      n++;
    }
    variance = sqrt(variance/(float)(n));
    fprintf(stderr,"OC %d: mean %f variance %f\n",oc,mean,variance);
  }
}

pushevent(struct data_str *ev,ev_stack **evs) 
//struct data_str  *ev;
//ev_stack        **evs;
{
  ev_stack *current_stack_pointer;

  if ((current_stack_pointer=(ev_stack*)malloc(sizeof(ev_stack)))==NULL) 
    usage("can't allocate eventstack pointer.");
  current_stack_pointer->prev_ev=*evs;
  *evs=current_stack_pointer;
  memcpy(&(current_stack_pointer->event),ev,datastr_size);
}

pulldown_events(ev_stack **evs)
//ev_stack **evs;
{
  ev_stack *current_stack_pointer;

  while (*evs != NULL) {
    current_stack_pointer=(*evs)->prev_ev;
    free(*evs);
    *evs=current_stack_pointer;
  }
}

usage(char *complaint)  
//char *complaint;
{
  fprintf(stderr,"\n%s\n\n",complaint);
  fprintf(stderr,"%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
"usage:medpict",
"      medpict [-b(burstmode)][-o <output biasfile>][-d <directory>] ",
"              [-f (lbox|astd|berlin)] ",
"              [-h(histmode)][-e(eventsearch)][-c(OCcorrection)]",
"              [-B <inputbias>] <filename1> <filename2> <filename3> .. ",
"",
"description: medpict reads in a number of FITS files and can perform any",
"             combination of the following:",
"             1) median bias calculation (and output)",
"             2) overclock variation correction",
"             3) event search",
"             4) bias subtracted histogram output",
"             5) input bias subtraction, (any of the above) and output of",
"                bias differences");

  fprintf(stderr,"%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
"MEDPICT works on a fairly simple principle. Given a list of input FITS",
"files and perhaps a specified bias frame, MEDPICT can perform most simple",
"operations. Most of these simple operations are useful for basic data",
"analysis with the X-Ray CCDs. If only a list of input FITS files are",
"named, then these are read in and a median frame is calculated. Optionally,",
"a bias frame may be specified to replace the median frame calculation.",
"Furthermore, if both input bias frame and output bias frames are specified,",
"then a median bias is calculated, and the difference between this frame and",
"the input bias frame is output. This can be useful if one is interested in",
"secular changes in the bias level over the entire device. ",
"Given that a bias frame is subtracted from each frame, MEDPICT may output",
"RV format event lists and/or imaging region pixel histograms, for each FITS",
"file named. Individual files with standard suffices, `evlist' & `pixhist'",
"are saved if a directory path is specified via the `-d' switch. ",
"If a <directory> is not specified, then the analysis products (EITHER evlist",
"OR pixel histograms - NOT BOTH) are placed on the standard output.",
"Incidentally, both standard formats (ASTRO-D and LBOX) are supported, and",
"may be specified on the command line using the `-f' flag. (Format support",
"means internal knowledge of `good' overclock sampling regions and extent of",
"imaging regions.)");
    exit(1);
}



void printerror( int status)
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
  
  char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
  if (status)
    fprintf(stderr, "\n*** Error occurred during program execution ***\n");
  
  fits_get_errstatus(status, status_str);   /* get the error description */
  fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);
  
  /* get first message; null if stack is empty */
  if ( fits_read_errmsg(errmsg) )
    {
      fprintf(stderr, "\nError message stack:\n");
      fprintf(stderr, " %s\n", errmsg);
      
      while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
        fprintf(stderr, " %s\n", errmsg);
    }
  
  exit( status );       /* terminate the program, returning error status */
}      
