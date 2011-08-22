/**/
extern double c_dXhyper( double *theta, int Mx[], int Ny[], int N1x[] );

/**/
extern double c_ddXhyper( double *theta, int Mx[], int Ny[], int N1x[] );

/**/
extern double c_dddXhyper( double *theta, int Mx[], int Ny[], int N1x[] );

/**/
extern double c_eXhyper( double *theta, int Mx[], int Ny[] );

/**/
extern double c_vXhyper( double *theta, int Mx[], int Ny[] );

/**/
extern void c_pmfXhyper( double *theta, int Mx[], int Ny[], double value[] );

/**/
extern int c_modeXhyper( double *theta, int Mx[], int Ny[] );
