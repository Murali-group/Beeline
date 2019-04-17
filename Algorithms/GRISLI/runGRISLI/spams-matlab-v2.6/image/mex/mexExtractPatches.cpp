/*!
 * \file
 *                toolbox imageRestor
 *
 *                by Julien Mairal
 *                julien.mairal@inria.fr
 *
 *                File mexExtractPatches.h
 * \brief mex-file, function mexExtractPatches
 * Usage: X = mexExtractPatches(I,n,step); */



#include <mexutils.h>
#include <image.h>

template <typename T>
inline void callFunction(mxArray* plhs[],
      const mxArray*prhs[], const int nrhs) {
   if (!mexCheckType<T>(prhs[0])) 
      mexErrMsgTxt("type of argument 1 is not consistent");
   if (mxIsSparse(prhs[0])) 
      mexErrMsgTxt("argument 1 should be full");

   T* prI = reinterpret_cast<T*>(mxGetPr(prhs[0]));
   const mwSize* dims=mxGetDimensions(prhs[0]);
   int nx=static_cast<int>(dims[0]);
   int ny=static_cast<int>(dims[1]);
   int dimsV=mxGetNumberOfDimensions(prhs[0]);
   int V=dimsV == 2 ? 1: (int)dims[2];

   int n=static_cast<int>(mxGetScalar(prhs[1]));
   int st = nrhs == 2 ? 1 : static_cast<int>(mxGetScalar(prhs[2]));

   int numPatchesX = st == 1 ? nx-n+1 : (nx-n+1)/st + 1;
   int numPatchesY = st == 1 ? ny-n+1 : (ny-n+1)/st + 1;
   INTM numPatches=numPatchesX*numPatchesY;
   plhs[0]=createMatrix<T>(n*n*V,numPatches);
   T* prX = reinterpret_cast<T*>(mxGetPr(plhs[0]));
   Matrix<T> X(prX,n*n*V,numPatches);
   Image<T> I(prI,ny,nx,V);

   I.extractPatches(X,n,st);

}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
   if (nrhs != 2 && nrhs != 3)
      mexErrMsgTxt("Bad number of inputs arguments");

   if (nlhs != 1) 
      mexErrMsgTxt("Bad number of output arguments");

   if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
      callFunction<double>(plhs,prhs,nrhs);
   } else {
      callFunction<float>(plhs,prhs,nrhs);
   }
}




