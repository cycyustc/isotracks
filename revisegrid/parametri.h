#ifndef PARAMETRI_H
#define PARAMETRI_H
//#include <stdio.h>
//#include <stdlib.h>
/************* new global parameters **********/
int itsacompfile = 0; // flag for dbert_comp
#define MAXSTRLEN 160 /* maximum length for scratch strings */
/****** random numbers: do not use ran0 (doesn't work) ******/
#define FRAND ran1   /* random number generation routine, see NR, 
			ran1 is probably best choice (small numbers) */
//int idum ; /* random number sequence, change to long if not ran1 */
/* #define double float */

/******* labels para estagios evolutivos ***********/
#define LABEL_DEFAULT -1
#define LABEL_PMS 0  
#define LABEL_MS 1  
#define LABEL_SUBGIANT 2  
#define LABEL_RGB 3  
#define LABEL_CHEB 4
#define LABEL_CHEBLM 4
#define LABEL_CHEBA 4
#define LABEL_CHEBB 5
#define LABEL_CHEBC 6
#define LABEL_EAGB 7
#define LABEL_TPAGB 8
#define LABEL_POSTAGB 9
#define LABEL_WD 10

/****** dichiarazioni ******/
#define	N_MET	1 //0 //15 //10 for parcol and then merge //15 for parsec	/* numero maximo di metalicidades */
#define	MET_NTRA	112	/* For new massive tracks numero massimo di tracce per metallicita' */
#define	MET_NGRU	10 //20	/* numero massimo di gruppi per metallicita' */
#define	TRA_NPUN	10000 //7720 //10000 for parcol//7720 for parsec	/* numero massimo di punti in traccia */
#define	TRA_NCRI	30	/* numero massimo di punti critici */
typedef	struct traccia
{
  int	npun, ncri, nm_zvar;
  float	mass_ini, mass_fin, mass_hecf,
    age_ini, age_fin;
  //quantitities along the track:
  float	age[TRA_NPUN]; 
  float	logte[TRA_NPUN], logl[TRA_NPUN];
  float mass_act[TRA_NPUN], mloss[TRA_NPUN], xcomp[TRA_NPUN][5]; //end
  char stage[TRA_NPUN][11];
  int label[TRA_NPUN];
  int	cri[TRA_NCRI], cri_zvar[TRA_NCRI]; //what is cri_zvar?
}	TRACCIA;	
typedef	struct	gruppo
{
  int	inf, sup;
}	GRUPPO;
typedef	struct	griglia
{
  int	ntra, ngru, i_zvar; //what is i_zvar?
  float	z, y, x, mh;
  TRACCIA	tra[MET_NTRA];
  GRUPPO	gru[MET_NGRU];	
  char	file[80];
}	GRIGLIA;
#endif /* PARAMETRI_H */
