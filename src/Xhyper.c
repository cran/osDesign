/***********************************************************************************************************/
/* Compile with:                                                                                           */
/*                                                                                                         */
/*   R CMD SHLIB Xhyper.c Xhyper_functions.c /home/haneussj/Research/Dissertation/C_Library/Smath/Smath.c  */
/*                                                                                                         */
/***********************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "R.h" // R memory io
#include "Rmath.h" // R math functions
#include "Rinternals.h" // R internal functions
/*#include "/home/haneussj/Research/Dissertation/C_Library/Smath/Smath.h"*/
#include "Xhyper.h"


/********************/
/* 2-LEVEL EXPOSURE */
/********************/
/**/

/** PMF **/
/**/
void dXhyper( double *theta, int Mx[], int Ny[], int N1x[], double *value )
{
  *value = c_dXhyper( theta, Mx, Ny, N1x );
  return;
}

/** PMF -- gradient **/
/**/
void ddXhyper( double *theta, int Mx[], int Ny[], int N1x[], double *value )
{
  *value = c_ddXhyper( theta, Mx, Ny, N1x );
  return;
}

/** PMF -- second derivative **/
/**/
void dddXhyper( double *theta, int Mx[], int Ny[], int N1x[], double *value )
{
  *value = c_dddXhyper( theta, Mx, Ny, N1x );
  return;
}

/** Mode **/
/**/
void modeXhyper( double *theta, int Mx[], int Ny[], int *value )
{
  *value = c_modeXhyper( theta, Mx, Ny );
}


/** Expectation **/
/**/
void eXhyper( double *theta, int Mx[], int Ny[], double *value )
{
  *value = c_eXhyper( theta, Mx, Ny );
}


/** Variance **/
/**/
void vXhyper( double *theta, int Mx[], int Ny[], double *value ){
  *value = c_vXhyper( theta, Mx, Ny );
  return;
}


/* Entire PMF: used for generating random deviates */
/**/
void pmfXhyper( double *theta, int Mx[], int Ny[], double value[] )
{
  c_pmfXhyper( theta, Mx, Ny, value );
  return;
}
