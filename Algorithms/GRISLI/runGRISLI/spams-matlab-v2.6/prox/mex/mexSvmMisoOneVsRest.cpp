
/* Software SPAMS v2.1 - Copyright 2009-2011 Julien Mairal 
 *
 * This file is part of SPAMS.
 *
 * SPAMS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SPAMS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SPAMS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mex.h>
#include <mexutils.h>
#include <svm.h>

// w=mexSvmMiso(y,X,tablambda,param);

template <typename T>
inline void callFunction(mxArray* plhs[], const mxArray*prhs[],
      const int nlhs) {
   if (!mexCheckType<T>(prhs[0])) 
      mexErrMsgTxt("type of argument 1 is not consistent");
   if (mxIsSparse(prhs[0])) 
      mexErrMsgTxt("argument 1 should not be sparse");

   if (!mexCheckType<T>(prhs[1])) 
      mexErrMsgTxt("type of argument 2 is not consistent");

   if (!mxIsStruct(prhs[2])) 
      mexErrMsgTxt("argument 3 should be a struct");

   T* pry = reinterpret_cast<T*>(mxGetPr(prhs[0]));
   const mwSize* dimsy=mxGetDimensions(prhs[0]);
   INTM my=static_cast<INTM>(dimsy[0]);
   INTM ny=static_cast<INTM>(dimsy[1]);
   Vector<T> y(pry,my*ny);

   T* prX = reinterpret_cast<T*>(mxGetPr(prhs[1]));
   const mwSize* dimsX=mxGetDimensions(prhs[1]);
   INTM p=static_cast<INTM>(dimsX[0]);
   INTM n=static_cast<INTM>(dimsX[1]);
   Matrix<T> X(prX,p,n);

   const int nclasses=y.maxval()+1;
   plhs[0]=createMatrix<T>(p,nclasses);
   T* prw=reinterpret_cast<T*>(mxGetPr(plhs[0]));
   Matrix<T> W(prw,p,nclasses);

   const int max_it = getScalarStructDef<int>(prhs[2],"max_it",1000*n);
   const T eps = getScalarStructDef<T>(prhs[2],"eps",0.001);
   int threads = getScalarStructDef<T>(prhs[2],"threads",-1);
   const T lambda = getScalarStruct<T>(prhs[2],"lambda");
   const bool accelerated = getScalarStructDef<T>(prhs[2],"accelerated",false);
   if (threads == -1) {
      threads=1;
#ifdef _OPENMP
      threads =  MIN(MAX_THREADS,omp_get_num_procs());
#endif
   } 
   threads=init_omp(threads);
   miso_svm_onevsrest(y,X,W,lambda,eps,max_it,accelerated); 
}

   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
      if (nrhs != 3)
         mexErrMsgTxt("Bad number of inputs arguments");

      if (nlhs != 1) 
         mexErrMsgTxt("Bad number of output arguments");

      if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
         callFunction<double>(plhs,prhs,nlhs);
      } else {
         callFunction<float>(plhs,prhs,nlhs);
      }
   }




