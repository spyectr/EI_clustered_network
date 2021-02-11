/*=================================================================
 * IntAuxPhi - computes integral inside transfer function for LIF neuron
 *
 * Input:   double ain, double bin
 * Output:  double out
 *
 * Compile inside Matlab cmd line as follows
 * 1) choose compiler:
 * >> mex -setup 
 * and pick gcc (option 1)
 * 2) compile:
 * >> mex IntAuxPhi.c
 * 3) use as a regular Matlab function e.g.
 * >> out=IntAuxPhi(1.23,2.24)
 * 
 * NOTE: if the integral in C is infinite for two consecutive iterations, 
 *       return 1.e+300, since Matlab cannot handle infinity.
 *
 * Copyright 2008-2011 The MathWorks, Inc.
 * Created by Luca Mazzucato on 2/18/14. *
 *=================================================================*/
/* file extension: .mexmaci64 on current machine macbookair */
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/************************\
 *    ROUTINE trapzd    *
 \************************/

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double),double a,double b,int n)
{
    double x,tnm,sum,del;
    static double s;
    int it,j;
    
    
    if(n==1){
        return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } else {
        for(it=1,j=1;j<n-1;++j) it<<=1;
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for(sum=0.0,j=1;j<=it;++j,x+=del) sum+=FUNC(x);
        s=0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}

/*-------
// QSIMP
//-------*/
#define EPS 1.0e-6   
#define JMAX 20     

double qsimp(double (*func)(double),double a,double b)
{
    double trapzd(double (*func)(double),double a,double b,int n);
    int j;
    double s,st,ost,os;
    /*printf("\n");*/
    ost=os=-1.0e30;
    for(j=1;j<=JMAX;++j){
        st=trapzd(func,a,b,j);
        s=(4.0*st-ost)/3.0;
        if(j>5)
/*             if (isinf(st) && isinf(ost)){ */
            if (st>1.0e+300 && ost>1.0e+300){
                s=1.0e+300; /* in Matlab, largest allowed number is realmax~1e+308 */
                return s;
            }
            if (fabs(s-os)<EPS*fabs(os)||(s==0.0 && os==0.0)) return s;
        os=s;
        ost=st;
        /*printf("--- # qsimp step: %d\r",j); //
        //fflush(stdout); */
    }
/*    
//    printf("\n\nError: too many steps in routine qsimp.\n");
//    return -1;   // error case */
}

/**********************************************************
 *                                                        *
 *                         phi(x):                        *
 *                                                        *
 *                      sqrt(pi)*x^2*(1+erf(x))           *
 *                                                        *
 **********************************************************/


#define a1  -1.26551223
#define a2  1.00002368
#define a3  .37409196
#define a4  .09678418
#define a5  -.18628806
#define a6  .27886807
#define a7  -1.13520398
#define a8 1.48851587
#define a9 -.82215223
#define a10 .17087277

double phi(double z)
{
    double t,ef,at;
    double w;
    w = fabs(z);
    t = 1./(1. + 0.5 * w);
    at=a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10))))))));
    ef=t*exp(at);
    if(z>0.) {
        ef = 2.*exp(w*w)-ef;
    }
    ef*=1.772453851;
    return ef;
}


/* The computational routine performs the integral given (ain,bin) extremes and returns out */
void compRoutine(double *ain,double *bin,double *out,int Nain)
{
    int i;
    
    for (i=1;i<=Nain;i++) {
        out[i-1]=qsimp(phi, ain[i-1], bin[i-1]);
    }
}



void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *ain, *bin, *out; /* pointers to input & output matrices*/
    size_t Main,Nain,Mbin,Nbin; /* dimensions of inputs */
    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }

    
    /* dimensions of input matrices */
    Main   = mxGetM(prhs[0]);  
    Nain   = mxGetN(prhs[0]);
    Mbin   = mxGetM(prhs[1]);  
    Nbin   = mxGetN(prhs[1]);
    
    if (Main != Mbin) {
        mexErrMsgIdAndTxt("MATLAB:IntAuxPhi_vec:matchdims",
                "Dimensions of input vectors do not match.");
    }
    if (Nain != Nbin) {
        mexErrMsgIdAndTxt("MATLAB:IntAuxPhi_vec:matchdims",
                "Dimensions of input vectors do not match.");
    }

    if (Main != 1) {
    mexErrMsgIdAndTxt("MATLAB:IntAuxPhi_vec:matchdims",
            "Inputs must be row vectors.");  
    }

    
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(1, Nain, mxREAL);
  /* Assign pointers to input and output. */
  ain = mxGetPr(prhs[0]);
  bin = mxGetPr(prhs[1]);
  /* Assign a pointer to the output. */
  out = mxGetPr(plhs[0]);
    /* compute integral */
   compRoutine(ain,bin,out,Nain);   
}