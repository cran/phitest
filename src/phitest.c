#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> 
#include <R_ext/Utils.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXTSP 2500
#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define ITMAX 100
#define JMAX 100
#define EPS 3.0e-8

void nrerror(char error_text[])
     /*  Numerical Recipes standard error handler  */
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

/***********************************************************************/
long double zbrent(long double (func)(long double, long double, long double, 
      long double, int), long double par1, long double par2, long double es, 
      int n,long double x1, long double x2, long double tol)
{
  int iter;
  long double a=x1,b=x2,c=x2,d,e,min1,min2;
  long double fa=(func)(a,par1,par2,es,n);
  long double fb=(func)(b,par1,par2,es,n),fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
    return -1.0;
  }
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;                           
      fc=fa;                       
      e=d=b-a;
    }
    if (fabsl(fc) < fabsl(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabsl(b)+0.5*tol;     
    xm=0.5*(c-b);

    if (fabsl(xm) <= tol1 || fb == 0.0) { 
      return b;
    }
    if (fabsl(e) >= tol1 && fabsl(fa) > fabsl(fb)) {
      s=fb/fa;                 
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;      
      p=fabsl(p);
      min1=3.0*xm*q-fabsl(tol1*q);
      min2=fabsl(e*q);

      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;                      
	d=p/q;
      } else {
	d=xm;                     
	e=d;
      }
    }
    else {
      d=xm;                 
      e=d;
    }
    a=b;                   
    fa=fb;
    if (fabsl(d) > tol1)      
      b += d;
    else 
      b += SIGN(tol1,xm);
    fb=(*func)(b,par1,par2,es,n);
  }
  nrerror("Maximum number of iterations exceeded in zbrent");										    
  return 0.0;            
}


/***********************************************************************/
long double rtsafe(long double (fn)(long double, long double, long double, long double, int), long double (dfn)(long double, long double, long double, long double, int), long double z1, long double z2, long double z3, int z4, long double x1, long double x2, long double xacc)
     /* uses the Newton-Raphson method to find the root of a function known to lie in the interval [x1, x2]  */
     /* method is applied until accuracy is within +/- xacc  */
     /* funcd() is a user supplied function that returns both the function value and it's derivative at the point x  */
{
  void nrerror(char error_text[]);
  int j;
  long double df,dx,dxold,f,fh,fl;
  long double nuisance;
  long double temp,xh,xl,rts;
 
  fl = fn(x1,z1,z2,z3,z4);
  fh = fn(x2,z1,z2,z3,z4);

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    nrerror("Root must be bracketed in rtsafe");

  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  }
  else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
    
  f = fn(rts,z1,z2,z3,z4);
  df = dfn(rts,z1,z2,z3,z4);

  for (j=1;j<=JMAX;j++) {

    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0) || (fabs(2.0*f)>fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
/*  printf("\n option 1, rts=%Lf",rts); */

      if (xl==rts) return rts;
    }
    else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
/*       printf("\n option 2, rts=%Lf", rts); */
      if (temp==rts) return rts;
    }
   
    if (fabs(dx) < xacc) return rts;
    f = fn(rts,z1,z2,z3,z4);
    df = dfn(rts,z1,z2,z3,z4);

    if (f < 0.0)
      xl=rts;
    else xh=rts;

/*     printf("\n rts = %Lf, xl = %Lf, xh = %Lf",rts,xl,xh); */
  }
  nrerror("Maximum number of iterations exceeded in rtsafe");
  return 0.0;
}

/***************************************************/
long double *get_vector(n)
int n;
{
long double* v;

if( n<0)
  fprintf(stderr,"Request for %d long doubles.\n",n);
if( n<=0)return NULL;

v = (long double *)malloc((n)*sizeof(long double));
if( !v)fprintf(stderr,"Can't obtain %d long doubles.\n",n);
return v;
}

int *get_vector2(v,n)
int n;
long double* v;
{
if( n<0)
  fprintf(stderr,"Request for %d long doubles.\n",n);
if( n<=0)return NULL;

v = (long double *)malloc((unsigned) n*sizeof(long double));
if( !v)fprintf(stderr,"Can't obtain %d long doubles.\n",n);
return 0;
}

int    *get_ivector( n)
int n;
{
int *v;

if( n<0)
  fprintf(stderr,"Request for %d ints detected.\n",n);
if( n<=0)return NULL;
v = (int *)malloc((unsigned) n*sizeof(int));
if( !v)fprintf(stderr,"Can't obtain %d ints.\n",n);
return v;
}

free_vector( v,n)
long double *v;
int    n;
{
if(n)free( (char*) (v) );
}

free_ivector( v,n)
int    *v;
int    n;
{
if(n)free( (char*) v );
}


/*********************************************************************/
int N_start(long double s)
{
int nstart;
if (s >= 1.0)
     nstart = 2;
else
     nstart = 3;
return nstart;
}


/*********************************************************************/
long double Lambda_start(long double s, long double alpha)
{
long double lambda;
if (s > 1.0) 
     lambda = 1.0/(s*(s - 1.0)) * (powl(alpha/2.0, 1.0 - s) - 1.0);
else if (s == 1.0) 
     lambda = -logl(alpha/2.0);
else if (s == 0.0) 
     lambda = 0.5*( (1-sqrtl(1.0-alpha))*logl(1-sqrtl(1.0-alpha)) 
                        + (1+sqrtl(1.0-alpha))*logl(1+sqrtl(1.0-alpha)) ); 
else 
     lambda = 0.5/(s*(1.0-s))*(2.0 - powl(1.0 -sqrtl(1.0-alpha),1.0-s) 
                        - powl(1.0+sqrtl(1.0-alpha),1.0-s) ); 
return lambda;
}





/***************************************************************/
long double cover_prob(n, a, b, lambda)
int n;
long double  *a, *b;   /* a, b of length n+1 */
long double lambda;
{
  int i,k,m,j,L,K,M,N;
  long double sum,delta,P,D,p,r,en;

  long double t[2*n+20], s[2*n+20], x[2*n+20], H[2*n+20], H1[2*n+20], coeff[2*n+20];
  int up[2*n+20], low[2*n+20];

/*   printf("n = %4d\n", n); */
/*   printf("l = %9.4f\n", l); */

  L=n; /* Since I've already taken care of when we have 0s and 1s, I can use the whole range here */
  K=0;
  M=n;

  for (k=1;k<=L;k++) {
    t[k]=a[k];
/*   printf("x[%d] = %9.4f\n", k,t[k]); */
  }
  t[n+1]=1;
  t[0]=0;
/*   printf("x[%d] = %9.4f\n", n+1,t[n+1]); */
  for (k=0;k<=M-1;k++) {
    s[k]=b[k+1];
/*   printf("y[%d] = %9.4f\n", k,s[k]); */
  }
  s[M]=1;
  if (s[M]=1) M--;
  N=M+L-K+1;

  x[0]=0;low[0]=-1;up[0]=1;
  i=0;k=K+1;
  for (j=1;j<=N;j++)
    {
      if (t[k]<s[i]) {x[j]=t[k];low[j]=i-1;up[j]=k;k++;}
      else {x[j]=s[i];low[j]=i;up[j]=k;i++;}
    }
  x[N+1]=1;up[N+1]=L+1;low[N+1]=M;

  H[0]=1; for (j=1;j<=N+1;j++) H[j]=0;

  for (k=1;k<=N+1;k++)
    {
      delta=n*(x[k]-x[k-1]);
      D=exp(-delta);

      coeff[0]=1;
      for (j=1;j<=up[k]-low[k-1];j++) coeff[j]=coeff[j-1]*delta/j;
      for (m=low[k]+1;m<=up[k]-1;m++)
	{
	  sum=0; for (j=0;j<=m-low[k-1]-1;j++) sum +=H[m-j]*coeff[j];
	  H1[m]=D*sum;
	}

      for (m=low[k]+1;m<=up[k]-1;m++) H[m]=H1[m];
    }

  /*I. P{     < \eta_n <       } */

  p=1.00;
  i=1;
  j=n;

  r = sqrtl(2*3.1415926*n);
  en = expl(1.0)/n;

  while (i <= j)
    {
      if (p<r) p *= (j--)*en;
      else p *= (i++)*en;
    }
  P=H[n]*p;
/*   printf("P = %9.4f\n",P); */

  return P;

}



/************************************************************/
long double fnKs(long double p, long double lambda, long double F_n, long double s, int  n)
     /* returns the function value fn at the point x */
     /* x is what i'm solving for */
{
  long double temp;
  if (s == 0.0){
  	temp = p*logl(p/F_n)+(1.0-p)*logl((1.0-p)/(1.0-F_n))-lambda;
  }
  else if (s == 1.0) {
  	temp = F_n*logl(F_n/p)+(1.0-F_n)*logl((1.0-F_n)/(1.0-p))-lambda;
}
  else {
      temp = (1.0 - powl(F_n,s)/powl(p,s-1.0)
                       - powl(1.0-F_n,s)/powl(1.0-p,s-1.0))/(s*(1.0-s))-lambda;
  }
  return temp;
}


/************************************************************/
long double dfnKs(long double p, long double lambda, long double F_n, long double s, int  n)
     /* returns the derivative df at the point x */
     /* x is what i'm solving for */
{
  long double temp;
  if (s == 0.0){
    temp = logl(p/F_n)-logl((1.0-p)/(1.0-F_n)); /* need to check this! */
  }
  else if (s == 1.0) {
	temp = (1.0 - F_n)/(1.0 - p) - F_n/p;   /* need to check this! */
}
  else {
    temp = 1.0/s*( powl(1.0-F_n,s)/powl(1.0-p,s) - powl(F_n,s)/powl(p,s) );
  }
  return temp;
}


/*********************************************************************/
long double PofS_n(long double lambda, long double alpha,
                           long double nuisance, long double s, int n)
{
  int i,j;
  long double *a,*b;
  long double F_n;
  long double answer;
  long double LL,UL;
  long double temp;
  int K;

  a = get_vector(n+2);
  b = get_vector(n+2);

  if( !a || !b )
    fprintf(stderr, "Not enough memory to handle n=%d", n);

  /* working with the b[] */

  /* finding the i where a[i]=0 */
  /* may eventually want to change this condition to the b[] */

  if (s >= 1.0)
    K = 0;
  else if (s < 1.0 && s != 0.0)
    K = (int) 1.0*n - n*powl(1.0-lambda*s*(1.0-s),1.0/s);
  else
    K= (int) n*(1.0-expl(-1.0*lambda));

  /* setting b[n-i+1]=1 for these i  */

  for (i=(n-K+1); i<=n; i++) {
    b[i] = 1.0;
  }

  /* setting the remaining b[i], 2 <= i <= n-K+1 */

  for (i=2; i<(n-K+1); i++) {
    LL= (i-1.0)/n;
    UL= 1.0;
    F_n = (i-1.0)/n;
   
    if (s == 2.0){
      j=n+1-i;
      b[i] = 1.0 - ((j + n*lambda) - sqrtl(2.0*n*j*lambda + n*n*lambda*lambda - 2.0*j*j*lambda))/(n+2.0*n*lambda);
    }
    else if (s == -1.0) {
      j=n+1-i;
      b[i] = 1.0 - (j*(1.0+2*lambda) - sqrt(2*lambda*j*(n+4*n*lambda-4*lambda*lambda*j-4*lambda*j-j+4*n*lambda*lambda)))/(n*(1.0+2*lambda));
    }
    else if (s == 0.5) {
      j=n+1-i;
      temp = ((4.0-lambda)*sqrt(1.0*j/(1.0*n)) - sqrt((1.0-1.0*j/(1.0*n))*(8*lambda-lambda*lambda)))/4.0;
      b[i] = 1.0 - temp*temp;
    }
    else
      b[i]=rtsafe(fnKs,dfnKs, lambda, F_n, s, n, LL, UL, 1e-16);
  }
 
  /* setting b[1] */

  LL = 1e-16;  /* for some reason, I can't make this 0.0 without a fuss! */
  UL = 1.0;

  F_n = 0.0;



  if (s < 1.0)
    b[1] = b[2];
  else if (s == 2.0)
    b[1] = 1.0 - ((n + n*lambda) - sqrtl(2.0*n*n*lambda + n*n*lambda*lambda - 2.0*n*n*lambda))/(n+2.0*n*lambda);
  else if (s ==1.0)
    b[1] = 1.0 - expl(-lambda);
  else
    b[1]=rtsafe(fnKs,dfnKs, lambda, F_n, s, n, LL,UL,1e-16);

/*   printf("\n\n\n b[1]=%Lf \n\n\n",b[1]);   */
/* computes the a[] */

  for (i=1; i<=n; i++){
    a[i]=1-b[n-i+1];
  }
  answer=cover_prob(n,a,b, lambda)-(1.0-alpha);

  free_vector(a, n+2);
  free_vector(b, n+2);
  return answer;
}








/*********************************************************************/
long double Pval(long double lambda, double s, int n)
{
  int i,j;
  long double *a,*b;
  long double F_n;
  long double answer;
  long double LL,UL;
  long double temp;
  int K;

  a = get_vector(n+2);
  b = get_vector(n+2);

  if( !a || !b )
    fprintf(stderr, "Not enough memory to handle n=%d", n);

  /* working with the b[] */

  /* finding the i where a[i]=0 */
  /* may eventually want to change this condition to the b[] */

  if (s >= 1.0)
    K = 0;
  else if (s < 1.0 && s != 0.0)
    K = (int) 1.0*n - n*powl(1.0-lambda*s*(1.0-s),1.0/s);
  else
    K= (int) n*(1.0-expl(-1.0*lambda));

  /* setting b[n-i+1]=1 for these i  */

  for (i=(n-K+1); i<=n; i++) {
    b[i] = 1.0;
  }

  /* setting the remaining b[i], 2 <= i <= n-K+1 */

  for (i=2; i<(n-K+1); i++) {
    LL= (i-1.0)/n;
    UL= 1.0;
    F_n = (i-1.0)/n;
   
    if (s == 2.0){
      j=n+1-i;
      b[i] = 1.0 - ((j + n*lambda) - sqrt(2.0*n*j*lambda + n*n*lambda*lambda - 2.0*j*j*lambda))/(n+2.0*n*lambda);
    }
    else if (s == -1.0) {
      j=n+1-i;
      b[i] = 1.0 - (j*(1.0+2*lambda) - sqrtl(2*lambda*j*(n+4*n*lambda-4*lambda*lambda*j-4*lambda*j-j+4*n*lambda*lambda)))/(n*(1.0+2*lambda));
    }
    else if (s == 0.5) {
      j=n+1-i;
      temp = ((4.0-lambda)*sqrtl(1.0*j/(1.0*n)) - sqrtl((1.0-1.0*j/(1.0*n))*(8*lambda-lambda*lambda)))/4.0;
      b[i] = 1.0 - temp*temp;
    }
    else
      b[i]=rtsafe(fnKs,dfnKs, lambda, F_n, s, n, LL, UL, 1e-16);
  }
 
  /* setting b[1] */

  LL = 1e-16;  /* for some reason, I can't make this 0.0 without a fuss! */
  UL = 1.0;

  F_n = 0.0;



  if (s < 1.0)
    b[1] = b[2];
  else if (s == 2.0)
    b[1] = 1.0 - ((n + n*lambda) - sqrtl(2.0*n*n*lambda + n*n*lambda*lambda - 2.0*n*n*lambda))/(n+2.0*n*lambda);
  else if (s ==1.0)
    b[1] = 1.0 - expl(-lambda);
  else
    b[1]=rtsafe(fnKs,dfnKs, lambda, F_n, s, n, LL,UL,1e-16);

/*   printf("\n\n\n b[1]=%Lf \n\n\n",b[1]);   */
/* computes the a[] */

  for (i=1; i<=n; i++){
    a[i]=1-b[n-i+1];
  }
  answer=1.0 - cover_prob(n,a,b, lambda);
/*   answer = s; */

  free_vector(a, n+2);
  free_vector(b, n+2);
  return answer;
}

/*********************************************************************/
void Pval_wrapper(double *lambda, double *s, int *n, double *scr)
{
  long double LAMBDA, S;
  int  N;
  LAMBDA = *lambda;
  S = *s;
  N = *n;

/*   *scr =  S; */
 *scr = Pval(LAMBDA, S, N);

}


/*********************************************************************/
void Bands(long double lambda, long double s, int n, double * scr)
{
  int i,j;
  long double *b;
  long double F_n;
  long double answer;
  long double LL,UL;
  long double temp;
  int K;

  b = get_vector(n+2);

  if( !b )
    fprintf(stderr, "Not enough memory to handle n=%d", n);

  /* working with the b[] */

  /* finding the i where a[i]=0 */
  /* may eventually want to change this condition to the b[] */

  if (s >= 1.0)
    K = 0;
  else if (s < 1.0 && s != 0.0)
    K = (int) 1.0*n - n*powl(1.0-lambda*s*(1.0-s),1.0/s);
  else
    K= (int) n*(1.0-expl(-1.0*lambda));

  /* setting b[n-i+1]=1 for these i  */

  for (i=(n-K+1); i<=n; i++) {
    b[i] = 1.0;
  }

  /* setting the remaining b[i], 2 <= i <= n-K+1 */

  for (i=2; i<(n-K+1); i++) {
    LL= (i-1.0)/n;
    UL= 1.0;
    F_n = (i-1.0)/n;
   
    if (s == 2.0){
      j=n+1-i;
      b[i] = 1.0 - ((j + n*lambda) - sqrtl(2.0*n*j*lambda + n*n*lambda*lambda - 2.0*j*j*lambda))/(n+2.0*n*lambda);
    }
    else if (s == -1.0) {
      j=n+1-i;
      b[i] = 1.0 - (j*(1.0+2*lambda) - sqrtl(2*lambda*j*(n+4*n*lambda-4*lambda*lambda*j-4*lambda*j-j+4*n*lambda*lambda)))/(n*(1.0+2*lambda));
    }
    else if (s == 0.5) {
      j=n+1-i;
      temp = ((4.0-lambda)*sqrtl(1.0*j/(1.0*n)) - sqrtl((1.0-1.0*j/(1.0*n))*(8*lambda-lambda*lambda)))/4.0;
      b[i] = 1.0 - temp*temp;
    }
    else
      b[i]=rtsafe(fnKs,dfnKs, lambda, F_n, s, n, LL, UL, 1e-16);
  }
 
  /* setting b[1] */

  LL = 1e-16;  /* for some reason, I can't make this 0.0 without a fuss! */
  UL = 1.0;

  F_n = 0.0;



  if (s < 1.0)
    b[1] = b[2];
  else if (s == 2.0)
    b[1] = 1.0 - ((n + n*lambda) - sqrtl(2.0*n*n*lambda + n*n*lambda*lambda - 2.0*n*n*lambda))/(n+2.0*n*lambda);
  else if (s ==1.0)
    b[1] = 1.0 - expl(-lambda);
  else
    b[1]=rtsafe(fnKs,dfnKs, lambda, F_n, s, n, LL,UL,1e-16);

 
  for (i=1; i<=n; i++){
    scr[i]= b[i];
  }

  free_vector(b, n+2);
  
}

/*********************************************************************/
void Bands_wrapper(double *lambda, double *s, int *n, double *scr)
{
  long double LAMBDA, S;
  int  N;
  LAMBDA = *lambda;
  S = *s;
  N = *n;

/*   *scr =  S; */
 Bands(LAMBDA, S, N, scr);

}


/*********************************************************************/
long double Quant(long double ALPHA, long double S, int N)
{

  int nvalues = N;
  long double *lambda;
  lambda=get_vector(nvalues+1);
  long double s = S;
  long double alpha = ALPHA;
  
  int n; 
  int nstart;
  
  nstart = N_start(s);
  lambda[nstart-1] = Lambda_start(s,alpha);
  
  /* If n <= 100, proceed with the normal program */
  
  if (nvalues <= 100) {

    for(n=nstart; n <= nvalues; n++) {

      long double answer, answer2;
      long double LL, UL;

      UL = lambda[n-1];
      LL = 1e-16;
      if (s <=1 && n <=15)
	UL = lambda[n-1]+1.0;
  
      answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);

      if (answer == -1.0) {
	UL = UL+lambda[n-1];
	answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      }
      lambda[n] = answer;
    }  
  }

  /* If 100 < n <= 500 */

  else if (nvalues <= 500) {

    for(n=nstart; n <= 100; n++) {

      long double answer, answer2;
      long double LL, UL;
      
      UL = lambda[n-1];
      LL = 1e-16;
      if (s <=1 && n <=15)
	UL = lambda[n-1]+1.0;
      
      answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      
      if (answer == -1.0) {
	UL = UL+lambda[n-1];
	answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      }
      
      lambda[n] = answer;
    }  
    
    long double answer, answer2;
    long double LL, UL;
    
    UL = lambda[100];
    LL = 1e-16;
    
    answer=zbrent(PofS_n,alpha,19.0,s,nvalues, LL, UL, 1e-26);
    
    if (answer == -1.0) {
      UL = UL+lambda[100];
      answer=zbrent(PofS_n,alpha,19.0,s,nvalues, LL, UL, 1e-26);
    }
    
    lambda[nvalues] = answer;
    
  }

  /* If 500 < n <= 1000 */

  else if (nvalues <= 1000) {

    for(n=nstart; n <= 100; n++) {

      long double answer, answer2;
      long double LL, UL;

      UL = lambda[n-1];
      LL = 1e-16;
      if (s <=1 && n <=15)
	UL = lambda[n-1]+1.0;
  
      answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      
      if (answer == -1.0) {
	UL = UL+lambda[n-1];
	answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      }
  
      lambda[n] = answer;
    }  
    
    long double answer, answer2;
    long double LL, UL;
    
    UL = lambda[100];
    LL = 1e-16;
    
    answer=zbrent(PofS_n,alpha,19.0,s,500, LL, UL, 1e-26);
    
    if (answer == -1.0) {
      UL = UL+lambda[100];
      answer=zbrent(PofS_n,alpha,19.0,s,500, LL, UL, 1e-26);
    }
    
    lambda[500] = answer;

    UL = lambda[500];
    LL = 1e-16;
    
    answer=zbrent(PofS_n,alpha,19.0,s,nvalues, LL, UL, 1e-26);
    
    if (answer == -1.0) {
      UL = UL+lambda[500];
      answer=zbrent(PofS_n,alpha,19.0,s,nvalues, LL, UL, 1e-26);
    }
    
    lambda[nvalues] = answer;

  }


 /* If n > 1000 */

  else {

    for(n=nstart; n <= 100; n++) {

      long double answer, answer2;
      long double LL, UL;
      
      UL = lambda[n-1];
      LL = 1e-16;
      if (s <=1 && n <=15)
	UL = lambda[n-1]+1.0;
      
      answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      
      if (answer == -1.0) {
	UL = UL+lambda[n-1];
	answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
      }
      
      lambda[n] = answer;
    }  
    
    long double answer, answer2;
    long double LL, UL;
    
    UL = lambda[100];
    LL = 1e-16;
    
    answer=zbrent(PofS_n,alpha,19.0,s,500, LL, UL, 1e-26);
    
    if (answer == -1.0) {
      UL = UL+lambda[100];
      answer=zbrent(PofS_n,alpha,19.0,s,500, LL, UL, 1e-26);
    }
    
    lambda[500] = answer;
    
    int stop;
    int i;
    stop = (int) nvalues/1000;
    
    for(i=1; i <= stop; i++) { 

	UL = lambda[500];
	if (i > 1)
	  UL = lambda[1000*(i-1)];
      LL = 1e-16;
      
      int temp;
      temp = 1000*i;

      answer=zbrent(PofS_n,alpha,19.0,s,temp, LL, UL, 1e-26);
      
      if (answer == -1.0) {
	UL = UL+UL;
	answer=zbrent(PofS_n,alpha,19.0,s,temp, LL, UL, 1e-26);
      }
      
      lambda[1000*i] = answer;
      
    }

    UL = lambda[1000*stop];
    LL = 1e-16;
    
    answer=zbrent(PofS_n,alpha,19.0,s,nvalues, LL, UL, 1e-26);
    
    if (answer == -1.0) {
      UL = UL+lambda[500];
      answer=zbrent(PofS_n,alpha,19.0,s,nvalues, LL, UL, 1e-26);
    }
    
    lambda[nvalues] = answer;  
  }
  

  long double temp;
  temp = lambda[nvalues]; 
  return temp;
  
  free_vector(lambda,nvalues+1);
  
}

/*********************************************************************/
long double Quant_orig(long double ALPHA, long double S, int N)
{

  int nvalues = N;
  long double *lambda;
  lambda=get_vector(nvalues+1);
  long double s = S;
  long double alpha = ALPHA;
 
  int n; 
  int nstart;

 /*  nstart = N_start(s); */
/*   lambda[nstart-1] = Lambda_start(s,alpha); */

  nstart=1000;
  lambda[nstart-1] = Lambda_start(s,alpha);

  for(n=nstart; n <= nvalues; n++) {

  long double answer, answer2;
  long double LL, UL;

  UL = lambda[n-1];
  LL = 1e-16;
if (s <=1 && n <=15)
  UL = lambda[n-1]+1.0;
  
  answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);

  if (answer == -1.0) {
    UL = UL+lambda[n-1];
    answer=zbrent(PofS_n,alpha,19.0,s,n, LL, UL, 1e-26);
  }
  
  lambda[n] = answer;
  }  
  long double temp;
  temp = lambda[nvalues]; 
  return temp;

  free_vector(lambda,nvalues+1);
  
}

/*********************************************************************/
void Quant_wrapper(double *alpha, double *s, int *n, double *scr)
{
  long double ALPHA, S;
  int  N;
  ALPHA = *alpha;
  S = *s;
  N = *n;

/*   *scr =  S; */
 *scr = Quant(ALPHA, S, N);

}
