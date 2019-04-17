/*!
 * \file
 *                toolbox imageRestor
 *
 *                by Julien Mairal
 *                julien.mairal@inria.fr
 *
 *                File mexExtractPatches.h
 * \brief mex-file, function mexExtractPatches
 * Usage: I2 = mexCombinePatches(X,I0,n,st,lambda,normalize);
 * lambda (0.0) and st (1) are optional  */



#include <mexutils.h>
#include <image.h>

template <typename T>
inline void callFunction(mxArray* plhs[],
      const mxArray*prhs[], const int nrhs) {
   if (!mexCheckType<T>(prhs[0])) 
      mexErrMsgTxt("type of argument 1 is not consistent");
   if (mxIsSparse(prhs[0])) 
      mexErrMsgTxt("argument 1 should be full");
   if (!mexCheckType<T>(prhs[1])) 
      mexErrMsgTxt("type of argument 2 is not consistent");
   if (mxIsSparse(prhs[1])) 
      mexErrMsgTxt("argument 2 should be full");

   T* prX = reinterpret_cast<T*>(mxGetPr(prhs[0]));
   const mwSize* dims=mxGetDimensions(prhs[0]);
   int sizePatch=static_cast<int>(dims[0]);
   int numPatches=static_cast<int>(dims[1]);
   
   T* prI = reinterpret_cast<T*>(mxGetPr(prhs[1]));
   const mwSize* dimsI=mxGetDimensions(prhs[1]);
   int nx=static_cast<int>(dimsI[0]);
   int ny=static_cast<int>(dimsI[1]);
   int dimsV=mxGetNumberOfDimensions(prhs[1]);
   int V=dimsV == 2 ? 1: (int)dimsI[2];

   int n=static_cast<int>(mxGetScalar(prhs[2]));
   T lambda = nrhs < 4 ? T() : static_cast<T>(mxGetScalar(prhs[3]));
   int st = nrhs < 5 ? 1 : static_cast<int>(mxGetScalar(prhs[4]));
   bool normalize = nrhs < 6 ? true : static_cast<bool>(mxGetScalar(prhs[5]));

   int numPatchesX = st == 1 ? nx-n+1 : (nx-n+1)/st + 1;
   int numPatchesY = st == 1 ? ny-n+1 : (ny-n+1)/st + 1;
   INTM numPatchesT=numPatchesX*numPatchesY;
   if (numPatchesT != numPatches) 
      mexErrMsgTxt("Number of patches is not consistent");
   if (n*n*V != sizePatch)
      mexErrMsgTxt("Size of patches is not consistent");

   Matrix<T> X(prX,sizePatch,numPatches);
   Image<T> I(prI,ny,nx,V);
   plhs[0]=createImage<T>(nx,ny,V);
   T* prI2 = reinterpret_cast<T*>(mxGetPr(plhs[0]));
   Image<T> I2(prI2,ny,nx,V);
   I2.copy(I);
   I2.combinePatches(X,lambda,st,normalize);

}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
   if (nrhs < 3 && nrhs > 6)
      mexErrMsgTxt("Bad number of inputs arguments");

   if (nlhs != 1) 
      mexErrMsgTxt("Bad number of output arguments");

   if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
      callFunction<double>(plhs,prhs,nrhs);
   } else {
      callFunction<float>(plhs,prhs,nrhs);
   }
}




