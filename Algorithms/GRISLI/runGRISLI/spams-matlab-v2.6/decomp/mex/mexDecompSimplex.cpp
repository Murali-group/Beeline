/* Software SPAMS v2.5 - Copyright 2009-2014 Julien Mairal
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

/*!
 * \file
 *                toolbox decomp
 *
 *                by Yuansi Chen and Julien Mairal
 *                julien.mairal@inria.fr
 *
 *                File mexDecompSimplex.cpp
 * \brief mex-file, function mexDecompSimplex
 * Usage: [alpha] = mexDecompSimplex(X,Z);
 * output a dictionary Z
 */

#include <mexutils.h>
#include <arch.h>

template <typename T>
inline void callFunction(mxArray* plhs[], const mxArray*prhs[],const int nrhs,
      const int nlhs) {
      if (!mexCheckType<T>(prhs[0]))
        mexErrMsgTxt("type of argument 1 is not consistent");
      if (mxIsSparse(prhs[0])) 
        mexErrMsgTxt("argument 1 should be full");
      if (!mexCheckType<T>(prhs[1]))
        mexErrMsgTxt("type of argument 1 is not consistent");
      if (mxIsSparse(prhs[1])) 
        mexErrMsgTxt("argument 1 should be full");

     T* prX = reinterpret_cast<T*>(mxGetPr(prhs[0]));
     const mwSize* dimsX=mxGetDimensions(prhs[0]);
     int m=static_cast<int>(dimsX[0]);
     int n=static_cast<int>(dimsX[1]);

     T* prZ = reinterpret_cast<T*>(mxGetPr(prhs[1]));
     const mwSize* dimsZ=mxGetDimensions(prhs[1]);
     int mZ=static_cast<int>(dimsZ[0]);
     int p=static_cast<int>(dimsZ[1]);
     bool computeXtX = getScalarStructDef<bool>(prhs[2],"computeXtX",true);
     int numThreads = getScalarStructDef<int>(prhs[2],"numThreads",-1);

     if (m != mZ) mexErrMsgTxt("argument sizes are not consistent"); 
     
     Matrix<T> X(prX,m,n);
     Matrix<T> Z(prZ,m,p);
     SpMatrix<T> alpha;
     
     decompSimplex<T>(X, Z, alpha, computeXtX,numThreads);
     convertSpMatrix(plhs[0],alpha.m(),alpha.n(),alpha.n(),
            alpha.nzmax(),alpha.v(),alpha.r(),alpha.pB());
}

   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
      if (nrhs != 3)
         mexErrMsgTxt("Bad number of inputs arguments");

      if (nlhs != 1)
         mexErrMsgTxt("Bad number of output arguments");

      if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
         callFunction<double>(plhs,prhs,nrhs,nlhs);
      } else {
         callFunction<float>(plhs,prhs,nrhs,nlhs);
      }
   }


