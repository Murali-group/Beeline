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
 *                toolbox dictLearn
 *
 *                by Yuansi Chen and Julien Mairal 
 *                julien.mairal@inria.fr
 *
 *                File mexArchetypalAnalysis.cpp
 * \brief mex-file, function mexArchetypalAnalysis
 * Usage: [Z ,[A],[B]] = mexArchetypalAnalysis(X,param,[Z0]);
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
   if (!mxIsStruct(prhs[1])) 
      mexErrMsgTxt("argument 2 should be struct");
   if (nrhs==3) {
      if (!mexCheckType<T>(prhs[2]))
         mexErrMsgTxt("type of argument 3 is not consistent");
      if (mxIsSparse(prhs[2])) 
         mexErrMsgTxt("argument 3 should be full");
   }

   Matrix<T> X;
   getMatrix(prhs[0],X);
   Matrix<T> Z0;
   int p;
   if (nrhs==3) {
      getMatrix(prhs[2],Z0);
      p=Z0.n();
   } else {
      p=getScalarStruct<int>(prhs[1],"p");
   }
   bool robust = getScalarStructDef<bool>(prhs[1],"robust",false);
   T epsilon = getScalarStructDef<T>(prhs[1],"epsilon",1e-3);
   bool computeXtX = getScalarStructDef<bool>(prhs[1],"computeXtX",true);
   int stepsFISTA = getScalarStructDef<int>(prhs[1],"stepsFISTA",3);
   int stepsAS = getScalarStructDef<int>(prhs[1],"stepsAS",50);
   int numThreads = getScalarStructDef<int>(prhs[1],"numThreads",-1);
   bool randominit = getScalarStructDef<bool>(prhs[1],"randominit",true);

   plhs[0]=createMatrix<T>(X.m(),p);
   Matrix<T> Z;
   getMatrix(plhs[0],Z);
   SpMatrix<T> A;
   SpMatrix<T> B;
   if (nrhs==3) {
      archetypalAnalysis<T>(X,Z0,Z, A,B,robust,epsilon,computeXtX,stepsFISTA,stepsAS, numThreads);
   } else {
      archetypalAnalysis<T>(X,Z,A,B,robust,epsilon,computeXtX,stepsFISTA,stepsAS,randominit,numThreads);
   }

   if (nlhs >= 2) 
      convertSpMatrix(plhs[1],A.m(),A.n(),A.n(),
            A.nzmax(),A.v(),A.r(),A.pB());
   if (nlhs ==3) 
      convertSpMatrix(plhs[2],B.m(),B.n(),B.n(),
            B.nzmax(),B.v(),B.r(),B.pB());
};

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
   if (nrhs != 2 && nrhs != 3)
      mexErrMsgTxt("Bad number of inputs arguments");

   if (nlhs != 1 && nlhs != 2 && nlhs != 3)
      mexErrMsgTxt("Bad number of output arguments");

   if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
      callFunction<double>(plhs,prhs,nrhs,nlhs);
   } else {
      callFunction<float>(plhs,prhs,nrhs,nlhs);
   }
}


