/*==========================================================
 * forceCalcC.c - Calculate Forces for Coulomb potential using adjacency matrix for scalling or thresholding
 *
 * inputs: positions (inMatrix)
 * adjacency matrix (inMatrix2)
 * charges (inMatrix3)
 * fundamental constants (k)
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
void distanceCalc(double *fx,double *fy,double *r,double *A, double *epsilon, double *sigma, mwSize n)
{
  double rij[2];
  double dr;
  double fmag;
  mwSize i;
  mwSize j;
  mwSize k;
  /* multiply each element y by x */
  for (i=0; i<n*n; i++) {
    //mexPrintf("%f %i \n", A[i], i);
    if (A[i] > 0) {
        j = i%n;
        k = i/n;
        //mexPrintf("%i %i \n", j, k);
        dr = sqrt(pow(r[j] - r[k],2) + pow(r[j+n]-r[k+n],2));
        rij[0] = (r[j] - r[k])/dr;
        rij[1] = (r[j+n]-r[k+n])/dr;       
        fmag = A[i]*48*epsilon[0]*(-1.0*pow(sigma[0],12)/pow(dr,13) + pow(sigma[0],6)/pow(dr,7)/2.0);
        //mexPrintf("%f %f %f %f %i %i\n", dr, fmag, rij[1], rij[2], j, k);
        fx[i] = fmag*rij[0];
        fy[i] = fmag*rij[1];
    }
    else{
        fx[i] = 0;
        fy[i] = 0;
    }

  }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *position;               /* Nx3 input matrix */
    double *adjacency;              /* NxN input matrix */
    //double *charge;              /* Nx1 input matrix charge */
    size_t ncols,nrows;                   /* size of matrix */
    double *fx;              /* NxN force components in the x direction */
    double *fy;             /* NxN force components in the y direction */
    double *epsilon;        /* Energy scale*/
    double *sigma;          /* Length scale */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Four inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required.");
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
    position = mxGetPr(prhs[0]);
    adjacency = mxGetPr(prhs[1]);
    epsilon = mxGetPr(prhs[2]);
    sigma = mxGetPr(prhs[3]);
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)nrows,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)nrows,mxREAL);
    //mexPrintf("Size of array in kilobytes: %.0f\n\n", inMatrix[0][0]);

    /* get a pointer to the real data in the output matrix */
    fx = mxGetPr(plhs[0]);
    fy = mxGetPr(plhs[1]);
    //outMatrix = inMatrix;
    /* call the computational routine */
    distanceCalc(fx,fy, position,adjacency,epsilon,sigma,(mwSize)nrows);
}
