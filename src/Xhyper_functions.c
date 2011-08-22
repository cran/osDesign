#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "R.h" // R memory io
#include "Rmath.h" // R math functions
#include "Rinternals.h" // R internal functions

#define c_max(a,b) (((a)>=(b))?(a):(b))
#define c_min(a,b) (((a)>=(b))?(b):(a))

/********************/
/* 2-LEVEL EXPOSURE */
/********************/
/**/

/** Recursive relation function; based on algorithm provided in Liao and Rosen (2001) **/
/**/
double c_r( double *theta, int Mx[], int Ny[], int index )
{
  double value = *theta * ((Mx[1] - index + 1) * (Ny[1] - index + 1)) / (index * (Mx[0] - Ny[1] + index) );
  return value;
}

double c_rr( double *theta, int Mx[], int Ny[], int index )
{
  double value = *theta * ((Mx[1] - index + 1) * (Ny[1] - index + 1)) / ((index-1) * (Mx[0] - Ny[1] + index) );
  return value;
}

double c_rrr( double *theta, int Mx[], int Ny[], int index )
{
  double value = *theta * ((Mx[1] - index + 1) * (Ny[1] - index + 1)) / ((index-2) * (Mx[0] - Ny[1] + index) );
  return value;
}

/** Mode **/
/**/
int c_modeXhyper( double *theta, int Mx[], int Ny[] )
{
  int lower, upper, value;
  double a, b, c, q_temp, root;

  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );

  a = *theta - 1;
  b = - (*theta * (Mx[1] + Ny[1] + 2)) - Mx[0] + Ny[1];
  c = *theta * (Mx[1] + 1) * (Ny[1] + 1);

  if( b < 0 ) q_temp =  -0.5 * (b - sqrt(pow(b,2) - 4*a*c));
  else q_temp =  -0.5 * (b + sqrt(pow(b,2) - 4*a*c));

  root = c / q_temp;
  value = (int) root;

  if( value < lower || value > upper )
    {
      root = q_temp / a;
      value = (int) root;
   }

  return value;
}


/** Probability Mass Function **/
/**/
double c_dXhyper( double *theta, int Mx[], int Ny[], int N1x[] ){
  int mode, lower, upper, i;
  double denom, numer, fi, value;

  mode = c_modeXhyper( theta, Mx, Ny );
  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );
  
  numer = 1;
  denom = 1;

  fi = 1;
  for( i = (mode + 1); i <= upper; i++ ){
      fi *= c_r(theta, Mx, Ny, i);
      denom += fi;
      if( i == N1x[1] ) numer = fi;
    }
  
  fi = 1;
  for( i = (mode - 1); i >= lower; i-- ){
      fi /= c_r(theta, Mx, Ny, i+1);
      denom += fi;
      if( i == N1x[1] ) numer = fi;
    }

  value = numer / denom;
  return value;
}


/** Probability Mass Function -- gradient **/
/**/
double c_ddXhyper( double *theta, int Mx[], int Ny[], int N1x[] ){
  int mode, lower, upper, i;
  double denom, numer, fi, value;

  mode = c_modeXhyper( theta, Mx, Ny );
  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );
  
  numer = 1;
  denom = 1;

  fi = 1;
  for( i = (mode + 1); i <= upper; i++ ){
      fi *= c_rr(theta, Mx, Ny, i);
      denom += fi;
      if( i == N1x[1] ) numer = fi;
    }
  
  fi = 1;
  for( i = (mode - 1); i >= lower; i-- ){
      fi /= c_rr(theta, Mx, Ny, i+1);
      denom += fi;
      if( i == N1x[1] ) numer = fi;
    }

  value = numer / denom;
  return value;
}



/** Probability Mass Function -- second derivative **/
/**/
double c_dddXhyper( double *theta, int Mx[], int Ny[], int N1x[] ){
  int mode, lower, upper, i;
  double denom, numer, fi, value;

  mode = c_modeXhyper( theta, Mx, Ny );
  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );
  
  numer = 1;
  denom = 1;

  fi = 1;
  for( i = (mode + 1); i <= upper; i++ ){
      fi *= c_rrr(theta, Mx, Ny, i);
      denom += fi;
      if( i == N1x[1] ) numer = fi;
    }
  
  fi = 1;
  for( i = (mode - 1); i >= lower; i-- ){
      fi /= c_rrr(theta, Mx, Ny, i+1);
      denom += fi;
      if( i == N1x[1] ) numer = fi;
    }

  value = numer / denom;
  return value;
}


/** Expectation **/
/**/
double c_eXhyper( double *theta, int Mx[], int Ny[] )
{
  int mode, lower, upper, i;
  double denom, numer, fi, value;
  
  mode = c_modeXhyper( theta, Mx, Ny );
  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );
  
  numer = mode;
  denom = 1;
  
  fi = 1;
  for( i = (mode + 1); i <= upper; i++ ){
      fi *= c_r(theta, Mx, Ny, i);
      numer += (i * fi);
      denom += fi;
    }
  
  fi = 1;
  for( i = (mode - 1); i >= lower; i-- ){
      fi /= c_r(theta, Mx, Ny, i+1);
      numer += (i * fi);
      denom += fi;
    }
  
  value = numer / denom;
  return value;
}


/** Variance **/
/**/
double c_vXhyper( double *theta, int Mx[], int Ny[] )
{
  int mode, lower, upper, i;
  double denom, numer_1, numer_2, fi, value;
  
  mode = c_modeXhyper( theta, Mx, Ny );
  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );
  
  numer_1 = mode;
  numer_2 = pow(mode, 2);
  denom = 1;

  fi = 1;
  for( i = (mode + 1); i <= upper; i++ ){
      fi *= c_r(theta, Mx, Ny, i);
      numer_1 += (i * fi);
      numer_2 += (pow(i, 2) * fi);
      denom += fi;
    }
  
  fi = 1;
  for( i = (mode - 1); i >= lower; i-- ){
      fi /= c_r(theta, Mx, Ny, i+1);
      numer_1 += (i * fi);
      numer_2 += (pow(i, 2) * fi);
      denom += fi;
    }
  
  value = (numer_2 / denom) - pow( numer_1 / denom, 2 );
  return value;
}


/* Entire PMF: used for generating random deviates */
/**/
void c_pmfXhyper( double *theta, int Mx[], int Ny[], double value[] )
{
  int mode, lower, upper, i;
  double denom, fi;

  mode = c_modeXhyper( theta, Mx, Ny );
  lower = c_max( 0, Ny[1] - Mx[0] );
  upper = c_min( Ny[1], Mx[1] );
  
  value[mode] = 1;
  denom = 1;

  fi = 1;
  for( i = (mode + 1); i <= upper; i++ ){
      fi *= c_r(theta, Mx, Ny, i);
      value[i] = fi;
      denom += fi;
    }
  
  fi = 1;
  for( i = (mode - 1); i >= lower; i-- ){
      fi /= c_r(theta, Mx, Ny, i+1);
      value[i] = fi;
      denom += fi;
    }

  for( i = (lower+1); i <= upper; i++ ){
    value[i] = value[i] / denom;
  }

  return;
}

