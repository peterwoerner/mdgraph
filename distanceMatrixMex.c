/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <stdio.h>


/* The computational routine */
void distanceCalc( double *z, double *y, double*x, mwSize n)
{
  mwSize i;
  mwSize j;
  mwSize k;
  /* multiply each element y by x */
  //mexPrintf('%d', n);
  for (i=0; i<n*n; i++) {
    j = i%n;
    k = i/n;
    
    z[i] = sqrt(pow(y[j] - x[k],2) + pow(y[j+n]-x[k+n],2));// + pow(y[j+2*n] - x[k+2*n],2));
    //mexPrintf("%i %i %f %f %f %f %f\n", j, k, y[j], x[k], y[j+n], x[k+n], z[i]);
  }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *inMatrix;               /* Nx3 input matrix */
    double *inMatrix2;              /* Nx3 input matrix */
    size_t ncols,nrows;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","One input required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input multiplier must be type double.");
    }
    
    /* make sure the second input argument is type double */
    /*if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    */
    /* check that number of rows in second input argument is 1 */
    /*if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }*/
    
    /* get the value of the scalar input  */
    //multiplier = mxGetScalar(prhs[0]);

    /* create a pointer to the real data in the input matrix  */
    inMatrix = mxGetPr(prhs[0]);
    inMatrix2 = mxGetPr(prhs[1]);
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)nrows,mxREAL);
    //mexPrintf("Size of array in kilobytes: %.0f\n\n", inMatrix[0][0]);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
   
    //outMatrix = inMatrix;
    /* call the computational routine */
    distanceCalc(outMatrix,inMatrix,inMatrix2,(mwSize)nrows);
}
