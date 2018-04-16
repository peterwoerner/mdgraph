/*==========================================================
 * forceCalcCList.c - Calculate Forces for Coulomb potential using adjacency list for scalling or thresholding
 * This should 
 * inputs: positions (inMatrix)
 * adjacency matrix (inMatrix2) list form of matrix with rows of index1, index2, scale
 * charges (inMatrix3)
 * fundamental constants (k)
 * 
 *
 * 
 *
 * This is a MEX-file for MATLAB.
 * 
 *
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <stdio.h>


/* The computational routine */
void distanceCalc(double *fx,double *fy,double *r,double *A,double *q, double *pf, mwSize n, mwSize length)
{
  double rij[2];
  double dr;
  double fmag;
  int i;
  mwSize j;
  mwSize k;
  /* multiply each element y by x */
  for (i=0; i<length; i++) {
    //mexPrintf("%f %i \n", A[i], i);       
    j = A[i];
    k = A[i+length];
    //mexPrintf("%i %i %i\n", i, j, k);
    dr = sqrt(pow(r[j] - r[k],2) + pow(r[j+n]-r[k+n],2));
    //mexPrintf("dr = %f\n",dr);
    rij[0] = (r[j] - r[k])/dr;
    rij[1] = (r[j+n]-r[k+n])/dr;  
    //mexPrintf("rij  = [%f,%f]\n",rij[0], rij[1]);
    //mexPrintf("q[j], q[k], A, k = [%f,%f,%f,%f]", q[j],q[k],A[i],pf[0]);
    fmag = -1*A[i+2*length]*pf[0]*q[j]*q[k]/dr/dr;
    //mexPrintf("%f %f %f %f %f %i %i\n", dr, fmag, rij[1], rij[2], j, k);
    fx[n*k+j] = fmag*rij[0];
    fx[n*j + k] = -1*fmag*rij[0]
    fy[n*k+j] = fmag*rij[1];
    fy[n*j + k] = -1*fmag*rij[1]

  }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *position;               /* Nx2 input matrix */
    double *adjacency;              /* mx3 input matrix */
    double *charge;              /* Nx1 input matrix charge */
    size_t ncols,nrows, length;                   /* size of matrix */
    double *fx;              /* NxN force components in the x direction */
    double *fy;             /* NxN force components in the y direction */
    double *k;  /* fundamental constant for Coulomb's laws (contains information about units). */

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
    charge = mxGetPr(prhs[2]);
    k = mxGetPr(prhs[3]);
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    length = mxGetM(prhs[1]);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)nrows,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)nrows,mxREAL);
    //mexPrintf("Size of array in kilobytes: %.0f\n\n", inMatrix[0][0]);

    /* get a pointer to the real data in the output matrix */
    fx = mxGetPr(plhs[0]);
    fy = mxGetPr(plhs[1]);
    //outMatrix = inMatrix;
    /* call the computational routine */
    distanceCalc(fx,fy, position,adjacency,charge,k,(mwSize)nrows, (mwSize)length);
}
