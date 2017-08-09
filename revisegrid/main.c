
/****************************************************/
/*	PREPARA.C	*/
/*	Prepara le tracce per l'interpolazione  */
/****************************************************/
//#define AGE_RESOLUTION 1.0e-7 //incluso per evitare due eta uguali
#define AGEMIN_RESOLUTION 40 // minimo numero di eta equidistanti
#define LOGL_RESOLUTION 0.04 //default resolution
#define LOGTE_RESOLUTION 0.01 //default resolution
#define INFITFACTOR4CBUR 2 // infitting factor for central burning stages
#define VERBOSE

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	"parametri.h"

/* #define VERBOSE */
/* #define VERBOSE_DEB */

GRIGLIA	low, hb, mas; //low[N_MET], hb[N_MET], mas[N_MET];
int	n_met, n_met_zams, n_met_mas, n_met_agb;
//float   eta_reimers ;
FILE	*fp_nomi, *fp_nomi2, *fp_log;
float logl_resolution = LOGL_RESOLUTION, 
  logte_resolution = LOGTE_RESOLUTION,
  agemin_resolution = AGEMIN_RESOLUTION;

/**************************** MAIN **********************************/
int	main(int argc, char *argv[])
{
  void leggi_griglie_tutte(); //, writes_interpola() ;
  char inputfile[MAXSTRLEN], outputfile[MAXSTRLEN];

  /* parses command line for options */
  if (argc<3) 
    {
      printf("Error: too few arguments\n");
      printf("Usage: main <inputfile> <outputfile> <res_logL> <res_logTe>\n");
      printf("       main <inputfile> <outputfile> (for default resolution of 0.04 0.01)\n");
      printf("       (0.02 0.005 used in PARSEC paper, 0.04 0.01 used in web version)\n");
      exit(1);
    } 
  else
    {
      strcpy(inputfile,argv[1]) ;
      strcpy(outputfile,argv[2]) ;
    }
  if (argc>3) 
    {
      logl_resolution = atof(argv[3]);
      logte_resolution = atof(argv[4]);
    }

  fp_nomi = fopen(inputfile, "r") ;
  fp_nomi2 = fopen(outputfile, "w") ;
  fp_log = fopen("logfile", "w") ;
  //printf ("\nModel input = %s ,", inputfile) ; 
  /* preparazione */
  leggi_griglie_tutte(fp_nomi, fp_nomi2);

  fclose (fp_nomi) ;
  fclose (fp_nomi2) ;
  fclose (fp_log) ;
  return 0;
}


/************************* LEGGI TRACCE *****************************/
void	leggi_griglie_tutte(fp_nomi, fp_nomi2)
     FILE	*fp_nomi, *fp_nomi2 ;
{
  char	nome[240], nome_low[240], nome_hb[240], 
    nome_tp[240], pwd[240], lixo[500] ;
  register int	iz, im, imtocopy;
  void	leggi_griglia(), writes_interpola(); //, leggi_griglia_agb()
  FILE	*fp_inp, *fp_out; //, *fp_outtrack ;
  float hb_massmax, low_massmax, mas_massmin, tpagb_massmin;
  int i, nlow, nhb, ii, itmp;
  int trac_max_lines = 0;
  int addstrackatend, lastpointtoadd, lastcritoadd, nnotpagbtoadd, 
    firstpointtoadd, firstcritoadd, ntrackinttosplit, ntracktoskip;
  float agetosubtract;
  int checktp;

  #define Z_X_SOLAR 0.0207
  //  kind_tpagb=0;
  //#define VERBOSE
  fscanf(fp_nomi, "%d", &n_met);
  fprintf(fp_nomi2, "%d\n", n_met);
  n_met_mas = n_met ;
  for (iz=0; iz<n_met; iz++)
    {
      itmp=-1;
      nnotpagbtoadd=0;
      checktp = 0 ;
      fscanf(fp_nomi, "%f %f", &low.z, &low.y);
      fscanf(fp_nomi, "%s", nome_low);
      fgets(lixo, 390, fp_nomi) ;
      printf("LOW %s\n",nome_low);
      fgets(lixo, 390, fp_nomi) ;
      sscanf(lixo, "%s", nome_hb);
      printf("HB  %s\n",nome_hb);
      sscanf(lixo, "%*s %s", nome_tp);
      fgets(lixo, 390, fp_nomi) ;
      sscanf(lixo, "%s", nome);
      printf("INT %s\n",nome);
      if (sscanf(lixo, "%*s %s", nome_tp)==1) {
	printf("TP %s\n",nome);
	checktp = 1 ;}
      low.x = 1.0 - low.z - low.y ;
      low.mh = log((low.z/low.x)/Z_X_SOLAR)/log(10.0);
      hb.z = low.z ;
      hb.y = low.y ;
      hb.x = low.x ;
      hb.mh = low.mh ;
      mas.z = low.z ;
      mas.y = low.y ;
      mas.x = low.x ;
      mas.mh = low.mh ;
     /* per le low : */
#ifdef	VERBOSE
      printf("\nLeggendo tracce LOW :\t%s", nome_low);
#endif
      fp_inp = fopen(nome_low, "r");
      leggi_griglia (fp_inp, &low, 0, nome_tp) ;
      for (im=0; im<low.ntra; im++)
	fscanf(fp_inp, "%f", &low.tra[im].mass_hecf);
      //leggi_interfaccia (fp_inp, &low[iz]) ; 
      fclose(fp_inp);
      nlow = low.ntra-1;
      low_massmax = low.tra[nlow].mass_ini;
      printf("largest mass LOW is :\t%f\n", low_massmax);
      /* per le hb : */
#ifdef	VERBOSE
      printf("\nLeggendo tracce HB :\t%s", nome_hb);
#endif
      fp_inp = fopen(nome_hb, "r");
      leggi_griglia (fp_inp, &hb, checktp, nome_tp) ;
      for (im=0, tpagb_massmin = 0.0; im<hb.ntra; im++)
	{
	  fscanf(fp_inp, "%f", &hb.tra[im].mass_hecf);
	  // stores min mass to attach tpagb later:
	  if (tpagb_massmin==0.0 && hb.tra[im].mass_hecf>0.0)
	    tpagb_massmin = hb.tra[im].mass_ini;
	  // imposes that higher masses will have interpolated tpagb:
	  if (tpagb_massmin>0.0 && hb.tra[im].mass_hecf==0.0)
	    hb.tra[im].mass_hecf = 1.0;}
      //leggi_interfaccia (fp_inp, &hb[iz]) ; 
      fclose(fp_inp);
      nhb = hb.ntra-1;
      hb_massmax = hb.tra[nhb].mass_ini;
      printf("largest mass HB is :\t%f\n", hb_massmax);
      /* per le intemass : */
#ifdef	VERBOSE
      printf("\nLeggendo tracce INTEMASS :\t%s", nome);
#endif
      fp_inp = fopen(nome, "r");
      leggi_griglia (fp_inp, &mas, checktp, nome_tp) ;
      for (im=0; im<mas.ntra; im++)
	fscanf(fp_inp, "%f", &mas.tra[im].mass_hecf);
      //leggi_interfaccia (fp_inp, &mas[iz]) ; 
      fclose(fp_inp);
      //mas_massmin = mas[iz].tra[0].mass_ini;
      //printf("lowest mass INT is :\t%f", mas_massmin);
      mas_massmin = mas.tra[ mas.gru[2].inf-1 ].mass_ini ;
      //printf(", but real lowest INT is :\t%f\n", mas_massmin);
      printf("lowest mass INT is :\t%f\n", mas_massmin);
      
      /*********************************************************************/
      strcpy(pwd, getenv("PWD"));
      /*************************replicates LOW******************************/
      strcat(nome_low, "2");
      fp_out = fopen(nome_low, "w");
      printf("Printing new LOW tracks to %s \n", nome_low);
     
      //////////////////////////////////
      if (mas_massmin >= 0.0)
      //if (mas[iz].tra[0].mass_ini > 0.0)
	{
	  //adds IMS track at the end
	  addstrackatend=1;
	  // last point of INT track is the one corresponding to last 
	  // critical point of the last LOW (supposedly the TRGB) 
	  ntrackinttosplit = mas.gru[2].inf;
	  //ntrackinttosplit = 1;
	  lastcritoadd = low.tra[ low.ntra-1 ].ncri-1;
	  //lastcritoadd = 5;
	  lastpointtoadd = mas.tra[ntrackinttosplit].cri[ lastcritoadd ];
	  //lastpointtoadd = 100;
	  agetosubtract = mas.tra[ntrackinttosplit].age[ lastpointtoadd ];
	  
	  if (mas_massmin == low_massmax)
	    {
	      im=low.ntra-1;
	      printf(" will replace INT splitted %dth (first %d points of %.2f Msun) to %dth already in LOW\n",
		     ntrackinttosplit, lastpointtoadd, mas_massmin, im);
	    }
	  else
	    {
	      im=low.ntra;
	      low.ntra++;
	      low.gru[low.ngru-1].sup++;
	      printf(" will append INT splitted %dth (%.2f Msun) to %dth in LOW\n",
		     ntrackinttosplit, mas_massmin, im);
	    }
	  low.tra[im].mass_ini = mas.tra[ntrackinttosplit].mass_ini;
	  low.tra[im].mass_hecf = 0.0;
	  low.tra[im].age_fin = mas.tra[ntrackinttosplit].age[lastpointtoadd];
	  low.tra[im].npun = lastpointtoadd+1;
	  low.tra[im].ncri = lastcritoadd+1;
	  for (ii=0; ii<=lastcritoadd; ii++)
	    low.tra[im].cri[ii] = mas.tra[ntrackinttosplit].cri[ii];
	  if (ii>=3) // to eliminate NEAR_ZAMS!!!!
	    {
	      //low[iz].tra[im].cri[3] = low[iz].tra[im].cri[4]-1;
	      //low[iz].tra[im].cri[2] = low[iz].tra[im].cri[4]-2;
	    }
	  for (ii=0; ii<=lastpointtoadd; ii++)
	    {
	      strcpy(low.tra[im].stage[ii], ""); //correct the labels
	      low.tra[im].age[ii] = mas.tra[ntrackinttosplit].age[ii]; 
	      low.tra[im].logl[ii] = mas.tra[ntrackinttosplit].logl[ii]; 
	      low.tra[im].logte[ii] = mas.tra[ntrackinttosplit].logte[ii]; 
	      if (itsacompfile==1)
		{
		  low.tra[im].mass_act[ii] = mas.tra[ntrackinttosplit].mass_act[ii]; 
		  low.tra[im].mloss[ii] = mas.tra[ntrackinttosplit].mloss[ii]; 
		  for (i=0; i<5; i++)
		    low.tra[im].xcomp[ii][i] = mas.tra[ntrackinttosplit].xcomp[ii][i]; 
		}
	    }
	  //correct the lables
	  for (ii=0; ii<=lastcritoadd; ii++)
	    strcpy(low.tra[im].stage[low.tra[im].cri[ii]],
		   low.tra[im-1].stage[low.tra[im-1].cri[ii]]);
	}
      else
	{
	  //if there is no IMS, does nothing
	  addstrackatend=0;
	}
      //////////////////////////////////

      /* prints numbers of groups, eliminating previous initial block  */
      fprintf(fp_out, "%d\n", low.ngru);
      for (i=0; i<low.ngru; i++)
	{
	  //if (addstrackatend==1 && i==low[iz].ngru-1) 
	  //  low[iz].gru[i].sup++;
	  fprintf(fp_out, "%d %d\n", low.gru[i].inf+1, low.gru[i].sup+1);
	}
      /* prints numbers and masses of tracks */
      fprintf(fp_out, "%d\n", low.ntra);
      for (i=0; i<low.ntra; i++)
	  fprintf(fp_out, "%d ", low.tra[i].npun); 
      fprintf(fp_out, "\n");
      for (i=0; i<low.ntra; i++)
	  fprintf(fp_out, "%f ", low.tra[i].mass_ini); 
      fprintf(fp_out, "\n");
      /* prints all tracks */
      for (im=0; im<low.ntra; im++)
	{
	  for (i=0; i<low.tra[im].npun; i++)
	    {
	      fprintf(fp_out, "%.12e %.5f %.5f ", 
		      low.tra[im].age[i], 
		      low.tra[im].logl[i], 
		      low.tra[im].logte[i]);
	      if (itsacompfile==1)
		{
		  fprintf(fp_out, "%.5f %.3g ", 
			  low.tra[im].mass_act[i], low.tra[im].mloss[i]);
		  for (ii=0; ii<5; ii++)
		    fprintf(fp_out, "%.3e ", low.tra[im].xcomp[i][ii]);		    
		}
	      if (i==0)
		fprintf(fp_out, "M=%f ", low.tra[im].mass_ini);
	      fprintf(fp_out, "%s %d ", 
		      low.tra[im].stage[i], i);
	      fprintf(fp_out, "\n");
	    }
	  if (low.tra[im].npun>trac_max_lines) 
	    trac_max_lines = low.tra[im].npun;
	}
      /* prints eqpoints for all tracks */
      for (i=0; i<low.ntra; i++)
	fprintf(fp_out, "%d ", low.tra[i].ncri); 
      fprintf(fp_out, "\n");
      for (im=0; im<low.ntra; im++)
	{
	  for (i=0; i<low.tra[im].ncri; i++)
	    fprintf(fp_out, "%d ", low.tra[im].cri[i]+1);
	  fprintf(fp_out, "\t!M=%f\n", low.tra[im].mass_ini);
	}
      for (i=0; i<low.ntra; i++)
	fprintf(fp_out, "%f ", low.tra[i].mass_hecf); 
      fprintf(fp_out, "\n");
      /* closes */
      fclose(fp_out);
      /****************************replicates HB******************************/
      strcat(nome_hb, "2");
      fp_out = fopen(nome_hb, "w");
      printf("Printing new HB tracks to %s \n", nome_hb);
     
      //////////////replace/add last HB track////////////////////
      if (addstrackatend==1) //whenever there are IMS tracks available
	{
	  ntrackinttosplit = mas.gru[2].inf;
	  firstcritoadd = low.tra[ low.ntra-1 ].ncri-1;
	  firstpointtoadd = mas.tra[ntrackinttosplit].cri[ firstcritoadd ];
	  agetosubtract = mas.tra[ntrackinttosplit].age[ firstpointtoadd ];
	  
	  if (mas_massmin == hb_massmax)
	    {
	      im=hb.ntra-1;
	      printf(" will replace INT splitted %dth (after %d point of %.2f Msun) to %dth already in HB\n",
		     ntrackinttosplit, firstpointtoadd, mas_massmin, im);
	    }
	  else
	    {
	      im=hb.ntra;
	      hb.ntra++;
	      hb.gru[hb.ngru-1].sup++;
	      printf(" will append INT splitted %dth (%.2f Msun) to %dth in HB\n",
		     ntrackinttosplit, mas_massmin, im);
	    }
	  hb.tra[im].mass_ini = mas.tra[ntrackinttosplit].mass_ini;
	  hb.tra[im].mass_hecf = mas.tra[ntrackinttosplit].mass_hecf;
	  hb.tra[im].age_fin = mas.tra[ntrackinttosplit].age_fin -
	    agetosubtract;
	  hb.tra[im].npun = mas.tra[ntrackinttosplit].npun -
	    firstpointtoadd;
	  hb.tra[im].ncri = mas.tra[ntrackinttosplit].ncri - 
	    firstcritoadd + 1;
	  for (ii=0; ii<hb.tra[im].ncri; ii++)
	    hb.tra[im].cri[ii] = 
	      mas.tra[ntrackinttosplit].cri[ii+firstcritoadd] -
	      firstpointtoadd ;  
	  //correction for lower number of cri in INT wrt HB:
	  hb.tra[im].cri[ hb.tra[im].ncri-1 ] = 
	    hb.tra[im].cri[ hb.tra[im].ncri-2 ];
	  hb.tra[im].cri[ hb.tra[im].ncri-2 ] = 
	    hb.tra[im].cri[ hb.tra[im].ncri-3 ] + 1;
	  for (ii=0; ii<hb.tra[im].npun; ii++)
	    {
	      strcpy(hb.tra[im].stage[ii], ""); //correct the labels
	      hb.tra[im].age[ii] = 
		mas.tra[ntrackinttosplit].age[ii+firstpointtoadd] -
		agetosubtract; 
	      hb.tra[im].logl[ii] = 
		mas.tra[ntrackinttosplit].logl[ii+firstpointtoadd]; 
	      hb.tra[im].logte[ii] = 
		mas.tra[ntrackinttosplit].logte[ii+firstpointtoadd]; 
	      if (itsacompfile==1)
		{
		  hb.tra[im].mass_act[ii] = mas.tra[ntrackinttosplit].mass_act[ii+firstpointtoadd]; 
		  hb.tra[im].mloss[ii] = mas.tra[ntrackinttosplit].mloss[ii+firstpointtoadd]; 
		  for (i=0; i<5; i++)
		    hb.tra[im].xcomp[ii][i] = mas.tra[ntrackinttosplit].xcomp[ii+firstpointtoadd][i]; 
		}
	    }
	  //correct the labels
	  for (ii=0; ii<hb.tra[im].ncri; ii++)
	    strcpy(hb.tra[im].stage[hb.tra[im].cri[ii]],
		   hb.tra[im-1].stage[hb.tra[im-1].cri[ii]]);
	}
      //////////////////////////////////

      ////////// duplicate/triplicate first track to have tpagb  ///////////
      for (im=0; im<hb.ntra; im++)
	if (hb.tra[im].mass_ini >= tpagb_massmin)
	  {
	    if (im==0) nnotpagbtoadd=2;
	    else nnotpagbtoadd=1;
	    itmp=im;
	    break;
	  }
      printf(" min mass for TPAGB is %.4f Msun, 1st track is %.4f Msun: %d will be added at %d\n",
	     tpagb_massmin, hb.tra[0].mass_ini, nnotpagbtoadd, itmp);
      // first moves all the next HB tracks up by nnotpagbtoadd,
      // then fills the hole, by replicating itmp nnotpagbtoadd times
      hb.ntra += nnotpagbtoadd;
      for (im=hb.ntra-1; im>=itmp; im--)
	{
	  if (im>=itmp+nnotpagbtoadd) 
	    imtocopy = im-nnotpagbtoadd;
	  else 
	    imtocopy = itmp;
	  hb.tra[im].npun = hb.tra[imtocopy].npun;
	  hb.tra[im].mass_ini = hb.tra[imtocopy].mass_ini;
	  for (i=0; i<hb.tra[im].npun; i++)
	    {
	      hb.tra[im].age[i] = hb.tra[imtocopy].age[i];
	      hb.tra[im].logl[i] = hb.tra[imtocopy].logl[i];
	      hb.tra[im].logte[i] = hb.tra[imtocopy].logte[i];
	      if (itsacompfile==1)
		{
		  hb.tra[im].mass_act[i] = hb.tra[imtocopy].mass_act[i];
		  hb.tra[im].mloss[i] = hb.tra[imtocopy].mloss[i];
		  for (ii=0; ii<5; ii++)
		    hb.tra[im].xcomp[i][ii] = hb.tra[imtocopy].xcomp[i][ii];
		}
	      strcpy( hb.tra[im].stage[i], hb.tra[imtocopy].stage[i]);
	    }
	  hb.tra[im].ncri = hb.tra[imtocopy].ncri;
	  for (i=0; i<hb.tra[im].ncri; i++)
	    hb.tra[im].cri[i] = hb.tra[imtocopy].cri[i];
	  //fixes the core mass:
	  if (im>=itmp+nnotpagbtoadd) 
	    hb.tra[im].mass_hecf = hb.tra[imtocopy].mass_hecf;
	  else 
	    hb.tra[im].mass_hecf = 0.0;
	  //fixes the inital mass to avoid deltaM=0.0 in fake interval:
	  if (nnotpagbtoadd==2)
	    hb.tra[0].mass_ini -= 0.000001;
	}
      // update groups: 
      hb.ngru++;
      for (i=hb.ngru-1; i>=1; i--)
	{
	  hb.gru[i].inf = hb.gru[i-1].inf+nnotpagbtoadd;
	  hb.gru[i].sup = hb.gru[i-1].sup+nnotpagbtoadd;
	}
      // limits change for 1st and 2nd interval:
      hb.gru[0].inf = 0;
      hb.gru[0].sup = itmp+nnotpagbtoadd-1;
      hb.gru[1].inf = itmp+nnotpagbtoadd; 
      //////////////////////////////////////////////////////////////////////

      /* prints numbers of groups, eliminating previous initial block  */
      fprintf(fp_out, "%d\n", hb.ngru);
      for (i=0; i<hb.ngru; i++)
	  fprintf(fp_out, "%d %d\n", hb.gru[i].inf+1, hb.gru[i].sup+1); 
      /* prints numbers and masses of tracks */
      fprintf(fp_out, "%d\n", hb.ntra);
      for (i=0; i<hb.ntra; i++)
	  fprintf(fp_out, "%d ", hb.tra[i].npun); 
      fprintf(fp_out, "\n");
      for (i=0; i<hb.ntra; i++)
	  fprintf(fp_out, "%f ", hb.tra[i].mass_ini); 
      fprintf(fp_out, "\n");
      /* prints all tracks */
      for (im=0; im<hb.ntra; im++)
	{
	  for (i=0; i<hb.tra[im].npun; i++)
	    {
	      fprintf(fp_out, "%.12e %.5f %.5f ", 
		      hb.tra[im].age[i], 
		      hb.tra[im].logl[i], 
		      hb.tra[im].logte[i]);
	      if (itsacompfile==1)
		{
		  fprintf(fp_out, "%.5f %.3g ", 
			  hb.tra[im].mass_act[i], hb.tra[im].mloss[i]);
		  for (ii=0; ii<5; ii++)
		    fprintf(fp_out, "%.3e ", hb.tra[im].xcomp[i][ii]);		    
		}
	      if (i==0)
		fprintf(fp_out, "M=%f ", hb.tra[im].mass_ini);
	      fprintf(fp_out, "%s %d ", 
		      hb.tra[im].stage[i], i);
	      fprintf(fp_out, "\n");
	    }
	  if (hb.tra[im].npun>trac_max_lines) 
	    trac_max_lines = hb.tra[im].npun;
	}
      /* prints eqpoints for all tracks */
      for (i=0; i<hb.ntra; i++)
	fprintf(fp_out, "%d ", hb.tra[i].ncri); 
      fprintf(fp_out, "\n");
      for (im=0; im<hb.ntra; im++)
	{
	  for (i=0; i<hb.tra[im].ncri; i++)
	    fprintf(fp_out, "%d ", hb.tra[im].cri[i]+1);
	  fprintf(fp_out, "\t!M=%f\n", hb.tra[im].mass_ini);
	}
      for (i=0; i<hb.ntra; i++)
	fprintf(fp_out, "%f ", hb.tra[i].mass_hecf); 
      fprintf(fp_out, "\n");
      /* closes */
      fclose(fp_out);

      /**************FIXES INT WHILE PRINTING****************************/
      strcat(nome, "2");
      // ******* STARTS PRINT INT FILE
      printf("Printing INT tracks, excluding low, to %s \n", nome);
      fp_out = fopen(nome, "w");
      //printf(" %d \n", mas[iz].ngru);
      //for (i=0; i<mas[iz].ngru; i++)
      //	  printf("%d %d\n", 
      //		  mas[iz].gru[i].inf-ntracktoskip, 
      //		  mas[iz].gru[i].sup-ntracktoskip); 
      ntracktoskip = mas.gru[1].sup+1;
      printf("Skipping %d tracks with < %.3f Msun from INT \n", 
	     ntracktoskip, mas.tra[ntracktoskip].mass_ini);
      /* prints numbers of groups, eliminating previous initial block  */
      fprintf(fp_out, "%d\n",  mas.ngru - low.ngru);
      for (i=0; i<mas.ngru; i++)
	if (i>=2)
	  fprintf(fp_out, "%d %d\n", 
		  mas.gru[i].inf-ntracktoskip+1, 
		  mas.gru[i].sup-ntracktoskip+1); 
      /* prints numbers of tracks, eliminating previous initial block  */
      fprintf(fp_out, "%d\n",  mas.ntra - low.ntra);
      //fprintf(fp_out, "%d ", tra.npun); 
      for (i=0; i<mas.ntra; i++)
	if (i>=ntracktoskip)
	  fprintf(fp_out, "%d ", mas.tra[i].npun); 
      fprintf(fp_out, "\n");
      //fprintf(fp_out, "%f ", tra.mass_ini); 
      for (i=0; i<mas.ntra; i++)
	if (i>=ntracktoskip)
	  fprintf(fp_out, "%f ", mas.tra[i].mass_ini); 
      fprintf(fp_out, "\n");
      /* prints composed track low+hb ************* OLD 
      for (i=0; i<tra.npun; i++)
	{
	  fprintf(fp_out, "%.12e %.5f %.5f ", 
		  tra.age[i], tra.logl[i], tra.logte[i]);
	  if (i==0)
	    fprintf(fp_out, "M=%f", tra.mass_ini);
	  fprintf(fp_out, "%s %d ", 
		  tra.stage[i], i);
	  fprintf(fp_out, "\n");
	}
      prints all other int tracks */
      for (im=1; im<mas.ntra; im++)
	{
	  if (im>=ntracktoskip)
	    for (i=0; i<mas.tra[im].npun; i++)
	      {
		fprintf(fp_out, "%.12e %.5f %.5f ", 
			mas.tra[im].age[i], 
			mas.tra[im].logl[i], 
			mas.tra[im].logte[i]);
		if (itsacompfile==1)
		  {
		    fprintf(fp_out, "%.5f %.3g ", 
			    mas.tra[im].mass_act[i], mas.tra[im].mloss[i]);
		    for (ii=0; ii<5; ii++)
		      fprintf(fp_out, "%.3e ", mas.tra[im].xcomp[i][ii]);		    
		  }
		if (i==0)
		  fprintf(fp_out, "M=%f ", mas.tra[im].mass_ini);
		fprintf(fp_out, "%s %d ", 
			mas.tra[im].stage[i], i);
		fprintf(fp_out, "\n");
	      }
	  if (mas.tra[im].npun>trac_max_lines) 
	    trac_max_lines = mas.tra[im].npun;
	}
      /* prints eqpoints for all tracks */
      //fprintf(fp_out, "%d ",  tra.ncri);
      for (i=0; i<mas.ntra; i++)
	if (i>=ntracktoskip)
	  fprintf(fp_out, "%d ", mas.tra[i].ncri); 
      fprintf(fp_out, "\n");

      //for (i=0; i<tra.ncri; i++)
      //fprintf(fp_out, "%d ", tra.cri[i]+1);
      //fprintf(fp_out, "\t!M=%f\n", tra.mass_ini);
      for (im=1; im<mas.ntra; im++)
	{
	  if (im>=ntracktoskip)
	    {
	      if (mas.tra[im].ncri>=3) // to eliminate NEAR_ZAMS!!!!
		{
		  //mas[iz].tra[im].cri[3] = mas[iz].tra[im].cri[4]-1;
		  //mas[iz].tra[im].cri[2] = mas[iz].tra[im].cri[4]-2;
		}
	      for (i=0; i<mas.tra[im].ncri; i++)
		fprintf(fp_out, "%d ", 
			mas.tra[im].cri[i]+1);
	      fprintf(fp_out, "\t!M=%f\n", mas.tra[im].mass_ini);
	    }
	}
      /* prints core masses */
      //fprintf(fp_out, "%f ", tra.mass_hecf); 
      for (i=0; i<mas.ntra; i++)
	if (i>=ntracktoskip)
	  fprintf(fp_out, "%f ", mas.tra[i].mass_hecf); 
      fprintf(fp_out, "\n");
 
      /* closes */
      fclose(fp_out);

      //print to list file
      strcpy(pwd, "isotrack/parsec");
      // ******* PRINTS PIECE OF FILELIST
      fprintf(fp_nomi2, "%f %f %f\n%s/%s\n%s/%s\n%s/%s\n", 
	      low.z, low.y, low.mh,
	      pwd, nome_low, pwd, nome_hb, pwd, nome);

    }
  printf("Max npun for resolution %f %f = %d\n",
	 logl_resolution, logte_resolution, trac_max_lines);
}




/**************** LEGGI GRIGLIA *******************/
void  leggi_griglia (fp_inp, grig, checktp, nome_tp)
     FILE	*fp_inp ;
     GRIGLIA	*grig ;
     int checktp ;
     char nome_tp[];
{
  void	leggi_traccia(), leggi_critici(), encurta_traccia();
  register int	im ;
  char	lixo[21];

  /* ECCO IL CASINO DEI PUNTI 0 0 !!!!!  */
  /* fscanf(fp_inp, "%*d %*d");  */
  fgets(lixo, 20, fp_inp);
  if (sscanf(lixo, "%*d %*d") == EOF) /* se c'e' un field solo */
    sscanf(lixo, "%d", &grig->ngru) ; /* leggi ngru da lixo e va avanti */
  else
    fscanf(fp_inp, "%d", &grig->ngru); /* leggi ngru dalla prossima riga */
#ifdef	VERBOSE
  fprintf(fp_log, "\n\t%d gruppi", grig->ngru);
#endif
  for (im=0; im<grig->ngru; im++)
    {
      fscanf(fp_inp, "%d %d", &grig->gru[im].inf, 
	     &grig->gru[im].sup);
      grig->gru[im].inf--, 
	grig->gru[im].sup--; 
    }
  fscanf(fp_inp, "%d", &grig->ntra);
#ifdef	VERBOSE
  fprintf(fp_log, "\n\t%d tracce", grig->ntra);
#endif
  for (im=0; im<grig->ntra; im++)
    fscanf(fp_inp, "%d", &grig->tra[im].npun);
  for (im=0; im<grig->ntra; im++)
    fscanf(fp_inp, "%f", &grig->tra[im].mass_ini);
  fgets(lixo, 20, fp_inp); //This seems necessary, because the first reading code in leggi_traccia is 'fgets'.
  for (im=0; im<grig->ntra; im++)
    {
      leggi_traccia( fp_inp, &grig->tra[im] ) ;
      grig->tra[im].age_fin = 
	grig->tra[im].age[ (grig->tra[im].npun-1) ] ;
    }
  for (im=0; im<grig->ntra; im++)
    {
      fscanf(fp_inp, "%d", &grig->tra[im].ncri);
#ifdef	VERBOSE
      fprintf(fp_log, "\n\t%d: %.2f Mo, %d punti, %d critici", im,
	     grig->tra[im].mass_ini, grig->tra[im].npun, 
	     grig->tra[im].ncri );
#endif
    }
  for (im=0; im<grig->ntra; im++)
    leggi_critici( fp_inp, &grig->tra[im] ) ;
  for (im=0; im<grig->ntra; im++)
    encurta_traccia(checktp, nome_tp, fp_inp, &grig->tra[im], &grig->z) ;
}

/************************ ENCURTA TRACCIA ********************/
void	encurta_traccia (checktp, nome_tp, fp_inp, trac, zeta)
     FILE	*fp_inp ;
     TRACCIA	*trac ;
     float zeta;
     int checktp ;
     char nome_tp[];
{
  register int	i, ii, ic, j;
  float loglres, logteres, ageres, s[TRA_NPUN], t[TRA_NPUN]; //ds, dt,
  char	lixo[180];
  int flag[TRA_NPUN];
  FILE *fp_tp;
  int ilast, itp, tpmodel, flagiflarger;
  float tpmass, tplogl, tplogte, tpmcore, tpage;

  if (checktp==1) // will open and store info for 1TP
    {
      tpmodel=0;
      //printf("mass %.3f : will check %d TPs in %s\n", 
      //     trac->mass_ini, checktp, nome_tp);
      fp_tp = fopen(nome_tp, "r");
      while (fscanf(fp_tp, "%d %*f %*f %*d %*f %f %f %f %*f %*f %f ", 
		   &itp, &tpmass, &tplogl, &tplogte, &tpmcore
		   ) != EOF) 
	{
	  for (ii=0; ii<25; ii++) fscanf(fp_tp, "%*f"); 
	  fscanf(fp_tp, "%f", &tpage); // that's the age at the correpsonding PARSEC track
	  fgets(lixo, 180, fp_tp);
	  //if (trac->mass_ini != tpmass) continue; 
	  if (trac->mass_ini < tpmass-1.0e-4 || trac->mass_ini > tpmass+1.0e-4) continue; 
	  tpmcore *= tpmass;
	  tpmodel = itp;
	  printf(" star %.4f found with mcore=%.4f, age=%.3e, model=%d\n",
		 tpmass, tpmcore, tpage, tpmodel);
	  break;
	}
      fclose(fp_tp);
    }
  ////return;
  s[0] = 0.0;
  t[0] = 0.0;
  ageres = (trac->age[trac->npun-1]-trac->age[0]) /
    AGEMIN_RESOLUTION;
  for (ii=1; ii<trac->npun; ii++)
    {
      loglres = logl_resolution;
      logteres = logte_resolution;
      // this doesn't work, need to redo:
      if (trac->label[ii]==LABEL_MS || 
	  trac->label[ii]==LABEL_CHEB ||
	  trac->label[ii]==LABEL_CHEBLM ||
	  trac->label[ii]==LABEL_CHEBA ||
	  trac->label[ii]==LABEL_CHEBB ||
	  trac->label[ii]==LABEL_CHEBC ) 
	{
	  loglres = logl_resolution / INFITFACTOR4CBUR;
	  logteres = logte_resolution / INFITFACTOR4CBUR;	  
	}
      s[ii] = s[ii-1] + 
	sqrt(
	     pow( (trac->logl[ii]-trac->logl[ii-1])/loglres, 2.0 ) + 
	     pow( (trac->logte[ii]-trac->logte[ii-1])/logteres, 2.0 )
	     ) ;
      t[ii] = t[ii-1] + (trac->age[ii]-trac->age[ii-1])/ageres ;
    }
  i=0;
  //ds=0.0;
  //dt=0.0;
  do // reject points too close to each other
    {
      flag[i] = 1 ;
      fprintf(fp_log, "aceita    ponto %d\n", i);
      for (ii=i+1; ii<trac->npun; ii++)
	if ( s[ii]-s[i]<1.0 && t[ii]-t[i]<1.0 )  
	  {
	    flag[ii] = 0 ;
	    fprintf(fp_log, "joga fora ponto %d < %d \n", ii, trac->npun-1);
	  }
	else
	  break;
      i=ii ;
    }
  while (i<trac->npun);
  for (ii=0; ii<trac->ncri; ii++) // re-accept critical points
    {
      flag[ trac->cri[ii] ] = 1 ;
      fprintf(fp_log, "aceita critico %d %d\n", trac->npun, trac->cri[ii]);
    }
  if (checktp==1 && tpmodel>0) // rejects points after 1TP, sets last point
    {
      ilast=trac->npun-1;
      flagiflarger=0;
      strcpy(lixo, trac->stage[ilast]) ; // captures name of last point in track
      for (ii=1; ii<trac->npun; ii++) // finds last possible age
	if (trac->age[ii] > tpage) // if the tracks lives longer than 1TP, will be cut
	  {
	    ilast=ii;
	    flagiflarger=1;
	    break; 
	  }
      if (flagiflarger==1) // will cut what remains of track
	for (ii=ilast+1; ii<trac->npun-1; ii++) // flags all rest as garbage
	  flag[ii] = 0 ;
      // now will fix last point of track, imposing it to be equal to 1TP
      flag[ilast-1] = 1 ; // will always include last useful point
      flag[ilast] = 1 ; // and also 1TP point
      trac->age[ilast] = tpage;
      trac->logl[ilast] = tplogl;
      trac->logte[ilast] = tplogte;
      //if (itsacompfile==1)
      //	{
      //	  trac->mass_act[ilast] = ;
      //	  trac->mloss[ilast] = ;
      //	  for (ii=0; ii<5; ii++)
      //	    trac->xcomp[ilast][ii] = ;		    
      //	}
      trac->mass_hecf = tpmcore;
      trac->npun = ilast+1;  
      trac->cri[ trac->ncri-1 ] = ilast; //set it as last cri one
      strcpy(trac->stage[trac->npun-1], lixo);
    }
  
  for (i=0, ii=0; i<trac->npun && i+ii<trac->npun; i++)
    {
      while (flag[i+ii]==0) ii++;
      trac->logl[i] = trac->logl[i+ii] ;
      trac->logte[i] = trac->logte[i+ii] ;
      trac->age[i] = trac->age[i+ii] ;
      if (itsacompfile==1)
	{
	  trac->mass_act[i] = trac->mass_act[i+ii] ;
	  trac->mloss[i] = trac->mloss[i+ii] ;
	  for (j=0; j<5; j++)
	    trac->xcomp[i][j] = trac->xcomp[i+ii][j] ;		    
	}
      strcpy(trac->stage[i], trac->stage[i+ii]) ;
      fprintf(fp_log, "%.3f Msun Z=%.4f: stored %d %d\n", trac->mass_ini, zeta,
	     i,i+ii);
      for (ic=0; ic<trac->ncri; ic++)
	if (trac->cri[ic]==i+ii) 
	  {
	    trac->cri[ic]=i ;
	    fprintf(fp_log, "critical %d %d\n",i,i+ii);
	  }
    }
  fprintf(fp_log, "initial %d ",trac->npun);
  trac->npun -= ii ;
  fprintf(fp_log, "stored %d rejected %d\n",trac->npun,ii); 
}

/************************ LEGGI TRACCIA ********************/
void	leggi_traccia (fp_inp, trac)
     FILE	*fp_inp ;
     TRACCIA	*trac ;
{
/* incluso per evitare due eta uguali: */
//#define AGE_RESOLUTION 1.0e-6

  register int	i, ii, ibug;
  char	lixo[380];
  for (ii=0; ii<trac->npun; ii++)
    {
      trac->stage[ii][0]=0; //null the stage to not be contaminated by the previous Z
      ibug=0;
      fgets(lixo, 380, fp_inp); //this will read first line, which can be a header
      // fscanf(fp_inp, "%f %f %f ", 
      //	     &trac->age[ii], &trac->logl[ii], &trac->logte[ii] );
      if (ii==0) // special procedure for first line
	{
	  sscanf(lixo, "%f %f %f", &trac->age[ii], &trac->logl[ii], &trac->logte[ii]);
	  //if (trac->logl[ii]==0.0)
	  //{
	  //ibug=1;
	  //fgets(lixo, 380, fp_inp); 
	  //}
	}
      if (strstr(lixo, "AGE") != NULL)  // in the case it's a header
	{
	  itsacompfile = 1;
	  fgets(lixo, 380, fp_inp); //reads next line
	}
      if (itsacompfile==0) //short track
	sscanf(lixo, "%f %f %f %*s %*s %s", 
	       &trac->age[ii], &trac->logl[ii], &trac->logte[ii],
	       trac->stage[ii]	     
	       ); 
      else
	sscanf(lixo, "%f %f %f %*s %*s %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %s", 
	       &trac->age[ii], &trac->logl[ii], &trac->logte[ii],
	       &trac->mass_act[ii],
	       &trac->xcomp[ii][0],
	       &trac->xcomp[ii][1],
	       &trac->xcomp[ii][2],
	       &trac->xcomp[ii][3],
	       &trac->xcomp[ii][4],
	       &trac->mloss[ii],
	       trac->stage[ii]	     
	       ); 
      //#ifdef	VERBOSE_DEB
      printf("\n\t\t%d %g %.4f %.4f", ii, 
	     trac->age[ii], trac->logl[ii], trac->logte[ii] );
      if (itsacompfile==1) //comp track
	{
	  for (i=0; i<5; i++)
	    trac->xcomp[ii][i] = pow(10.0, trac->xcomp[ii][i]);
	  printf(" %.4f %.4f", trac->mass_act[ii], trac->xcomp[ii][0] );
	}
      //#endif
      // to avoid bug in parsec files:
      if (trac->age[ii]==0.0 && ii>1) // to avoid bug in parsec files:
	{
	  //ii--;
	  ibug=1;
	}
      //WARNING: WAS THIS STILL NECESSARY?
      //if (ibug==1)  strcpy(trac->stage[ii],"MS_BEG");
    }
}

void	leggi_critici (fp_inp, trac)
     FILE	*fp_inp ;
     TRACCIA	*trac ;
{
  register int	ii;
  char	lixo[80];
  for (ii=0; ii<trac->ncri; ii++)
    {
      fscanf(fp_inp, "%d", &trac->cri[ii] );
      trac->cri[ii]-- ;
#ifdef	VERBOSE
      fprintf(fp_log, "\n\t\t%d %d", ii, trac->cri[ii] );
#endif
    }
  fgets(lixo, 80, fp_inp);

  if (trac->ncri==14) // to avoid bug in parsec INT files:
    {
      trac->ncri++;
      for (ii=trac->ncri-1; ii>4; ii--)
	trac->cri[ii]=trac->cri[ii-1];
      trac->cri[4]=trac->cri[3]+1;
    }
  if (trac->ncri==5) // to avoid bug in parsec LOW files:
    {
      trac->ncri++;
      for (ii=trac->ncri-1; ii>4; ii--)
	trac->cri[ii]=trac->cri[ii-1];
      trac->cri[4]=trac->cri[3]+1;
    }
}

//not used
void	leggi_critici_zvar (fp_inp, trac)
     FILE	*fp_inp ;
     TRACCIA	*trac ;
{
  register int	ii;
  char	lixo[80];
  for (ii=0; ii<trac->ncri; ii++)
    {
      fscanf(fp_inp, "%d", &trac->cri_zvar[ii] );
      trac->cri[ii]-- ;
#ifdef	VERBOSE
      fprintf(fp_log, "\n\t\t%d %d", ii, trac->cri_zvar[ii] );
#endif
    }
  fgets(lixo, 80, fp_inp);
}

//################################################
void writes_interpola(fp_out, gr1, gr2)
     FILE	*fp_out ;
     GRIGLIA	*gr1, *gr2 ;
{
  int i, ig, n1, n2, temp, temp1;
  
  for (ig=0; ig<gr1->ngru; ig++)
    {
      n1 = gr1->gru[ig].sup - gr1->gru[ig].inf;
      n2 = gr2->gru[ig].sup - gr2->gru[ig].inf;
      fprintf (fp_log, "interval %d: z1 from %d to %d, z2 from %d to %d \n",
	      ig, gr1->gru[ig].inf+1, gr1->gru[ig].sup+1, 
	      gr2->gru[ig].inf+1, gr2->gru[ig].sup+1 );
      for (i=0; i<=n1; i++) 
	if (n2==n1) /* discards nothing */
	  {
	    fprintf (fp_log, "n2==n1: %d in Z1=%.4f will be matched to %d in Z2=%.4f\n", 
		    gr1->gru[ig].inf+i+1, gr1->z, gr2->gru[ig].inf+i+1, gr2->z);
	    fprintf(fp_out, "%d ", gr2->gru[ig].inf + i + 1);
	  }
	else if (n2>n1) /* discards tracks in z2 when needed */
	  {
	    temp = (int)( 0.5+(float)i*(float)n2/(float)n1 ) - i ;
	    fprintf (fp_log, "n2>n1: %d in Z1=%.4f will be matched to %d in Z2=%.4f, %d\n", 
		    gr1->gru[ig].inf+i+1, gr1->z, gr2->gru[ig].inf+i+1+temp, gr2->z, 
		    temp);
	    fprintf (fp_out, "%d ", gr2->gru[ig].inf + i + 1 + temp );
	  }
	else if (n1>n2) /* discards tracks in z1 when needed */
	  {
	    if (i==0)
	      {
		fprintf (fp_out, "%d ", gr2->gru[ig].inf+1 );
		fprintf (fp_log, "n1>n2: %d in Z1 will be matched to %d in Z2\n", 
			gr1->gru[ig].inf+1, gr2->gru[ig].inf+1);
	      }
	    else
	      {
		temp1 = (int)( 0.5+(float)(i-1)*(float)n2/(float)n1 ) - (i-1) ;
		temp = (int)( 0.5+(float)i*(float)n2/(float)n1 ) - i ;
		if (temp!=temp1)
		  {
		    fprintf (fp_out, "0 " );
		    fprintf (fp_log, "n1>n2: %d in Z1 will not be matched in Z2, %d\n", 
			    gr1->gru[ig].inf+i+1, temp);
		  }
		else
		  {
		    fprintf (fp_out, "%d ", gr2->gru[ig].inf + i + 1 + temp );
		    fprintf (fp_log, "n1>n2: %d in Z1 will be matched to %d in Z2, %d\n", 
			    gr1->gru[ig].inf+i+1, gr2->gru[ig].inf+i+1+temp, 
			    +temp);
		  }
	      }
	  }
      fprintf(fp_log, "\n");
    }
  fprintf(fp_out, "\n");
}
