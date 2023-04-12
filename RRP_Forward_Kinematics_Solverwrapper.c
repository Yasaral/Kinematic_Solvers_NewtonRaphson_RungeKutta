

#include <math.h>  
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#define  _NAN      9.99999999999999E+305

void     eqns1     (double integtau, double VAR[], double VARp[], char boundary);
void     solve     (int size, double *a[], double b[], double x[]);
int      kutta     (void(*eqns)(double, double*, double*, char),
                    int numy, double y[], double *t, double integstp,
                    double abserr, double relerr, char com);


double   s1,s2,s3,theta,x,y;
double   Pi=3.141592653589793,z[94],_COEF[6][6],*COEF[6],RHS[6],RHSGUESS[6],GUESS[6];
double   r=151.28, TiltA_forward = 0, q1, q2, q3;
/* ................................ MAIN ............................. */
void   RRP_Forward_Kinematics_wrapper  (const double* posi, const double* posd_old, double* posd)
{
int      iloop;
double   integtau,taufinal=1.0,integstp = 0.25,abserr = 1.0e-08,relerr = 1.0e-07;
double   VAR[6];

/* Get inputs -- independent generalized coordinates */
q1 = posi[0];
q2 = posi[1];
q3 = posi[2];

/* Get guesses from the previous iteration */
x     = posd_old[0];
y     = posd_old[1];
theta = posd_old[2];
s1    = posd_old[3];
s2    = posd_old[4];
s3    = posd_old[5];


/* Initialize COEF pointers to beginning of each row */
for(iloop=0;  iloop<6;  iloop++)  COEF[iloop] = &(_COEF[iloop][0]);



/* Evaluate constants */
  z[1] = cos(TiltA_forward);
  z[2] = sin(TiltA_forward);
  z[3] = cos(q1);
  z[4] = sin(q1);
  z[5] = cos(q2);
  z[6] = sin(q2);
  z[7] = cos(q3);
  z[8] = sin(q3);
  z[59] = z[1]*z[3] - z[2]*z[4];
  z[60] = z[1]*z[4] + z[2]*z[3];
  z[65] = z[1]*z[5] - z[2]*z[6];
  z[66] = z[1]*z[6] + z[2]*z[5];
  z[75] = z[1]*z[7] - z[2]*z[8];
  z[76] = z[1]*z[8] + z[2]*z[7];


/* Initialize independent variable and integrator array */
  integtau = 0.0;
  VAR[0] = x;
  VAR[1] = y;
  VAR[2] = theta;
  VAR[3] = s1;
  VAR[4] = s2;
  VAR[5] = s3;

/* Store guess(es) and value(s) of function(s) at guess(es) */
kutta(eqns1,6,VAR,&integtau,integstp,abserr,relerr,0);
for(iloop=0;  iloop<6;  iloop++)
  { GUESS[iloop] = VAR[iloop];  RHSGUESS[iloop] = RHS[iloop]; }

/* Numerically integrate to obtain solution */
kutta(eqns1,6,VAR,&integtau,integstp,abserr,relerr,0);
while( integtau + .01*integstp < taufinal &&
       (iloop = kutta(eqns1,6,VAR,&integtau,integstp,abserr,relerr,1)) != 0);

	   
posd[0] = x;
posd[1] = y;
posd[2] = theta;	
posd[3] = s1;
posd[4] = s2;
posd[5] = s3;
   
	   
	   
	   
}


/* ................................ EQNS1 ............................. */
void     eqns1        (double integtau, double VAR[], double VARp[], char boundary)
{

/* Update variables with new values */
  x = VAR[0];
  y = VAR[1];
  theta = VAR[2];
  s1 = VAR[3];
  s2 = VAR[4];
  s3 = VAR[5];

  z[9] = cos(theta);
  z[10] = sin(theta);
  z[52] = 0.5000000000000003*z[9] - 0.8660254037844386*z[10];
  z[53] = 0.5000000000000003*z[10] + 0.8660254037844386*z[9];
  z[55] = z[1]*z[52] - z[2]*z[53];
  z[56] = z[1]*z[53] + z[2]*z[52];
  z[62] = z[2]*z[10] - z[1]*z[9];
  z[63] = -z[1]*z[10] - z[2]*z[9];
  z[64] = z[1]*z[10] + z[2]*z[9];
  z[68] = 0.5*z[9] + 0.8660254037844386*z[10];
  z[69] = 0.5*z[10] - 0.8660254037844386*z[9];
  z[71] = z[1]*z[68] - z[2]*z[69];
  z[72] = z[1]*z[69] + z[2]*z[68];
  z[78] = -0.5000000000000003*z[1]*(z[10]+1.732050807568876*z[9]) - 0.5000000000000003*
  z[2]*(z[9]-1.732050807568876*z[10]);
  z[79] = s1*z[78];
  z[80] = x + s1*z[55] - r*z[59];
  z[81] = 0.5000000000000003*z[1]*(z[9]-1.732050807568876*z[10]) - 0.5000000000000003*
  z[2]*(z[10]+1.732050807568876*z[9]);
  z[82] = s1*z[81];
  z[83] = y + s1*z[56] - r*z[60];
  z[84] = s2*z[64];
  z[85] = x + s2*z[62] - r*z[65];
  z[86] = s2*z[62];
  z[87] = y + s2*z[63] - r*z[66];
  z[88] = 0.5*z[1]*(1.732050807568877*z[9]-z[10]) - 0.5*z[2]*(z[9]+1.732050807568877*
  z[10]);
  z[89] = s3*z[88];
  z[90] = x + s3*z[71] - r*z[75];
  z[91] = 0.5*z[1]*(z[9]+1.732050807568877*z[10]) + 0.5*z[2]*(1.732050807568877*
  z[9]-z[10]);
  z[92] = s3*z[91];
  z[93] = y + s3*z[72] - r*z[76];

  COEF[0][0] = -1;
  COEF[0][1] = 0;
  COEF[0][2] = -z[79];
  COEF[0][3] = -z[55];
  COEF[0][4] = 0;
  COEF[0][5] = 0;
  COEF[1][0] = 0;
  COEF[1][1] = -1;
  COEF[1][2] = -z[82];
  COEF[1][3] = -z[56];
  COEF[1][4] = 0;
  COEF[1][5] = 0;
  COEF[2][0] = -1;
  COEF[2][1] = 0;
  COEF[2][2] = -z[84];
  COEF[2][3] = 0;
  COEF[2][4] = -z[62];
  COEF[2][5] = 0;
  COEF[3][0] = 0;
  COEF[3][1] = -1;
  COEF[3][2] = -z[86];
  COEF[3][3] = 0;
  COEF[3][4] = -z[63];
  COEF[3][5] = 0;
  COEF[4][0] = -1;
  COEF[4][1] = 0;
  COEF[4][2] = -z[89];
  COEF[4][3] = 0;
  COEF[4][4] = 0;
  COEF[4][5] = -z[71];
  COEF[5][0] = 0;
  COEF[5][1] = -1;
  COEF[5][2] = -z[92];
  COEF[5][3] = 0;
  COEF[5][4] = 0;
  COEF[5][5] = -z[72];
  RHS[0] = z[80];
  RHS[1] = z[83];
  RHS[2] = z[85];
  RHS[3] = z[87];
  RHS[4] = z[90];
  RHS[5] = z[93];
  solve(6,COEF,RHSGUESS,VARp);

}

int      kutta        ( void(*eqns)(double, double*, double*, char),
                        int numy, double y[], double *t,
                        double integstp, double abserr, double relerr, char com )
{
static double  f0[100], f1[100], f2[100], y1[100], y2[100];
static int     numcuts = 20;                     /* Max # cuts of integstp  */
static double  hc = 0;                           /* Last value of stepsize  */
char           entry = 1;                        /* Just entered routine    */
char           stepdouble;                       /* Double the stepsize     */
double         tfinal = *t + integstp;           /* Time at end of full step*/
double         tt, h;
int            i;

if(numy >= 100) {printf("\nERROR: INCREASE THE STATIC DOUBLE ARRAY SIZE IN kutta.\n" ); return 0;}
if( !com ) { (*eqns)(*t,y,f0,1);  return 1;}     /* Fill array f0 and return*/
if( numy == 0)  { hc = integstp;  return 1;}     /* Check for initial entry */
if( integstp == 0)  return 0;                    /* Cannot integrate forward*/
if( hc*integstp < 0 ) hc = -hc;                  /* Integrate backward      */
else if( hc == 0 )    hc = integstp;             /* Maybe initial entry     */
h  = hc;                                         /* Current stepsize        */
tt = *t + h;                                     /* Terminal time this step */
*t = tfinal;                                     /* Return updated t value  */

beginning:
while( tt+h != tt )                              /* Check round-off problems*/
  {
  double h2 = h * 0.5;                                     /* Half    of h  */
  double h3 = h / 3.0;                                     /* Third   of h  */
  double h6 = h / 6.0;                                     /* Sixth   of h  */
  double h8 = h * 0.125;                                   /* Eighth  of h  */
  if( com==2 || entry)
 {(*eqns)( tt-h,     y, f0, 1 );   entry=0; }              /* Boundary here */
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h3*f0[i];
  (*eqns)( tt-2*h3, y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h6*(f0[i] + f1[i]);
  (*eqns)( tt-2*h3, y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h8*(f0[i] + 3*f1[i]);
  (*eqns)( tt-h2,   y1, f2, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h2*(f0[i] - 3*f1[i] + 4*f2[i]);
  (*eqns)( tt,      y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y2[i] = y[i] + h6*(f0[i] + 4*f2[i] + f1[i]);
  stepdouble = 1;                                /* Assume need to double   */
  for(i=0;  i<numy;  i++)                        /* Check all equations     */
    {
    double error = fabs( y1[i] - y2[i] ) * 0.2;  /* Error in estimate       */
    double test  = fabs( y1[i] ) * relerr;       /* For relative error      */
    if( error >= test  &&  error >= abserr )     /* Error criterion not met?*/
      {                                          /* Restart w/ half stepsize*/
      hc = h = h2;                               /* Halve the stepsize      */
      tt -= h2;                                  /* Change terminal time    */
      if( numcuts-- > 0 )  goto beginning;       /* Back to beginning       */
      printf("\n THE STEPSIZE HAS BEEN HALVED TOO MANY TIMES; T = %15.8E"
             "\n NUMERICAL INTEGRATION FAILED TO CONVERGE.\n", *t=tt-h );
      (*eqns)(*t,y,f0,0);                        /* Fill for error display  */
      return 0;
      }
    if( stepdouble && 64*error>test && 64*error>abserr ) stepdouble = 0;
    }
  for(i=0;  i<numy;  i++)  y[i] = y2[i];         /* Update y for next step  */
  if( stepdouble && fabs(h+h)<=fabs(integstp) && fabs(tt+h+h)<=fabs(tfinal) )
    {hc=(h+=h);  numcuts++;}                     /* Double the stepsize     */
  if( tt == tfinal )                             /* End of integration      */
    { (*eqns)(tfinal,y,f0,2);  return 1;}        /* Derivatives at tfinal   */
  tt += h;                                       /* Update terminal time    */
  /*** If next jump puts tt past or very close to tfinal, adjust h and tt ***/
  if( (h>0 && tt>tfinal-0.1*h) || (h<0 && tt<tfinal-0.1*h) )
    { h = tfinal-(tt-h);  tt = tfinal; }         /* Adjust for last jump    */
  if( com == 1 )                                 /* Approx this derivative  */
    for(i=0;  i<numy;  i++)  f0[i] = f1[i];
  if( com == 3 )                                 /* Calc derivative once    */
    (*eqns)( tt-h, y, f0, 0 );
 }
printf("\nTHE STEPSIZE OF %15.7E IS TOO SMALL RELATIVE TO THE TERMINAL TIME"
       "\nOF %15.7E.  INTEGRATION HALTED BECAUSE OF NUMERICAL ROUND-OFF.\n"
       "\nTHE STEPSIZE MAY HAVE BEEN CUT TOO MANY TIMES.\n\n", h, *t=tt);
(*eqns)(*t,y,f0,0);                              /* Fill for error display  */
return 0;
}
 
void           solve  (int n, double *A[], double B[], double X[])
{
static double  scales[100];                      /* Row scaling factors     */
int            i, j, k, swapi;
double         rowmax, big, size, pivotReciprocal, ratio, sum, *swapa, swapx;
                                                 
if(n >= 100) {printf("\nERROR: INCREASE THE STATIC DOUBLE ARRAY SIZE IN solve.\n" ); exit(0);}
for(i=0;  i<n;  i++)                             /*** Begin decomposition ***/
  {
  for(rowmax=0.0, j=0;  j<n;  j++)               /* Check for zero row      */
    if( rowmax < fabs(A[i][j]) )  rowmax = fabs(A[i][j]);
  if( rowmax == 0.0 )                            /* Zero row found          */
    {printf("ALL ELEMENTS IN ROW %d OF COEF ARE ZEROS\n", i+1);    exit(0);}
  scales[i] = 1.0 / rowmax;                      
  X[i] = B[i];                                   /* Leave B matrix unchanged*/
  }
for(k=0;  k<n-1;  k++)                           /* Until 1 column is Left  */
  {
  for(big=0.0, i=k;  i<n;  i++)                  /* Check remaining rows    */
    {
    size = fabs( A[i][k] ) * scales[i];          /* Relative size this col  */
    if( size > big )  {swapi=i;   big=size;}     /* Largest relative entry  */
    }
  if( big == 0.0 )                               /* Zero pivot found        */
    {
    printf("A PIVOT ELEMENT ENCOUNTERED IN THE DECOMPOSITION OF COEF IS ZERO"
           "\nCOEFFICIENT MATRIX IS SINGULAR\n");
    exit(0);
    }
  if( swapi != k )                               /* Switch rows of A and X  */
    {
    swapa=A[k];  A[k]=A[swapi];  A[swapi]=swapa; /* Change row pointers     */
    swapx=X[k];  X[k]=X[swapi];  X[swapi]=swapx; /* Change X values         */
    }
  pivotReciprocal = 1.0 / A[k][k];               /* Value of pivot          */
  for(i=k+1;  i<n;  i++)                         /* Change lower rows       */
    {
    ratio = A[i][k] * pivotReciprocal;           /* Multiplicative factor   */
    X[i] -= ratio * X[k];                        /* Modify elements of X    */
    for(j=k+1;  j<n;  j++)  A[i][j] -= ratio * A[k][j];  
    }
  }
if( A[n-1][n-1] == 0.0 )                         /* Check last pivot        */
  {
  printf("A PIVOT ELEMENT ENCOUNTERED IN THE DECOMPOSITION OF COEF IS ZERO"
         "\nCOEFFICIENT MATRIX IS SINGULAR\n");
  exit(0);
  }
for(i=n-1;  i>=0;  i--)                          /* Solve upper * X = Z     */
  {
  for(sum=0.0, j=i+1;  j<n;  j++)  sum += A[i][j] * X[j];
  X[i] = (X[i] - sum) / A[i][i];
  }
}



