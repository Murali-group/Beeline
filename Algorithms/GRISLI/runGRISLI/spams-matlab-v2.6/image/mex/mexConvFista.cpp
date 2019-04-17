#include <mex.h>
#include <mexutils.h>
#include <fista.h>
#include <image.h>

// alpha = mexFistaFlat(I,D,alpha0,param)

using namespace FISTA;

template <typename T>
inline void callFunction(mxArray* plhs[], const mxArray*prhs[],
      const int nlhs) {
   if (!mexCheckType<T>(prhs[0])) 
      mexErrMsgTxt("type of argument 1 is not consistent");
   if (mxIsSparse(prhs[0])) 
      mexErrMsgTxt("argument 1 should not be sparse");

   if (!mexCheckType<T>(prhs[1])) 
      mexErrMsgTxt("type of argument 1 is not consistent");
   if (mxIsSparse(prhs[1])) 
      mexErrMsgTxt("argument 1 should not be sparse");

   if (!mexCheckType<T>(prhs[2])) 
      mexErrMsgTxt("type of argument 3 is not consistent");
   if (mxIsSparse(prhs[2])) 
      mexErrMsgTxt("argument 3 should not be sparse");

   if (!mxIsStruct(prhs[3])) 
      mexErrMsgTxt("argument 4 should be a struct");

   T* prI = reinterpret_cast<T*>(mxGetPr(prhs[0]));
   int V= mxGetNumberOfDimensions(prhs[0]);
   const mwSize* dimsX=mxGetDimensions(prhs[0]);
   const int h=static_cast<int>(dimsX[0]);
   const int w=static_cast<int>(dimsX[1]);
   const int nim= (V == 2 || static_cast<int>(dimsX[2])==1) ? 1 : static_cast<int>(dimsX[2]);
   Matrix<T> I(prI,h*w*nim,1);

   const mwSize* dimsD=mxGetDimensions(prhs[1]);
   int m=static_cast<int>(dimsD[0]);
   int p=static_cast<int>(dimsD[1]);
   T* prD = reinterpret_cast<T*>(mxGetPr(prhs[1]));
   Matrix<T> D(prD,m,p);
   ConvolveDictionary<T> convD(D,w,h,nim);

   int sizeEdge=static_cast<int>(sqrt(static_cast<T>(m/nim)));

   T* pr_alpha0 = reinterpret_cast<T*>(mxGetPr(prhs[2]));
   const mwSize* dimsAlpha=mxGetDimensions(prhs[2]);
   int sizeMapXY=(static_cast<int>(dimsAlpha[1]))*static_cast<int>(dimsAlpha[2]);
   int pAlpha=static_cast<int>(dimsAlpha[0]);
   Matrix<T> alpha0(pr_alpha0,sizeMapXY*pAlpha,1);

   plhs[0]=createImage<T>(p,(h-sizeEdge+1),(w-sizeEdge+1));
   //plhs[0]=createImage<T>(p,sizeMapXY,1);
   T* pr_alpha=reinterpret_cast<T*>(mxGetPr(plhs[0]));
   Matrix<T> alpha(pr_alpha,sizeMapXY*pAlpha,1);
   alpha.copy(alpha0);

   FISTA::ParamFISTA<T> param;
   param.num_threads = getScalarStructDef<int>(prhs[3],"numThreads",-1);
   param.max_it = getScalarStructDef<int>(prhs[3],"max_it",1000);
   param.tol = getScalarStructDef<T>(prhs[3],"tol",0.000001);
   param.it0 = getScalarStructDef<int>(prhs[3],"it0",100);
   param.pos = getScalarStructDef<bool>(prhs[3],"pos",false);
   param.linesearch_mode= getScalarStructDef<int>(prhs[3],"linesearch_mode",0);
   param.compute_gram = false;
   param.max_iter_backtracking = getScalarStructDef<int>(prhs[3],"max_iter_backtracking",1000);
   param.L0 = getScalarStructDef<T>(prhs[3],"L0",1.0);
   param.lambda= getScalarStructDef<T>(prhs[3],"lambda",1.0);
   param.lambda2= getScalarStructDef<T>(prhs[3],"lambda2",0.0);
   param.lambda3= getScalarStructDef<T>(prhs[3],"lambda3",0.0);
   param.admm = false; 
   param.lin_admm = false;
   getStringStruct(prhs[3],"regul",param.name_regul,param.length_names);
   param.regul = regul_from_string(param.name_regul);
   if (param.regul==INCORRECT_REG)
      mexErrMsgTxt("Unknown regularization");
   getStringStruct(prhs[3],"loss",param.name_loss,param.length_names);
   param.loss = loss_from_string(param.name_loss);
   if (param.loss==INCORRECT_LOSS)
      mexErrMsgTxt("Unknown loss");

   param.intercept = false;
   param.resetflow = getScalarStructDef<bool>(prhs[3],"resetflow",false);
   param.verbose = getScalarStructDef<bool>(prhs[3],"verbose",false);
   param.clever = getScalarStructDef<bool>(prhs[3],"clever",false);
   param.ista= getScalarStructDef<bool>(prhs[3],"ista",false);
   param.subgrad= getScalarStructDef<bool>(prhs[3],"subgrad",false);
   param.log= getScalarStructDef<bool>(prhs[3],"log",false);
   param.a= getScalarStructDef<T>(prhs[3],"a",T(1.0));
   param.b= getScalarStructDef<T>(prhs[3],"b",0);

   if (param.log) {
      mxArray *stringData = mxGetField(prhs[3],0,"logName");
      if (!stringData) 
         mexErrMsgTxt("Missing field logName");
      int stringLength = mxGetN(stringData)+1;
      param.logName= new char[stringLength];
      mxGetString(stringData,param.logName,stringLength);
   }

   if (param.num_threads == -1) {
      param.num_threads=1;
#ifdef _OPENMP
      param.num_threads =  MIN(MAX_THREADS,omp_get_num_procs());
#endif
   } 

   if (param.regul==GRAPH || param.regul==GRAPHMULT) 
      mexErrMsgTxt("Error: penalty not supported");
   if (param.regul==TREE_L0 || param.regul==TREEMULT || param.regul==TREE_L2 || param.regul==TREE_LINF) 
      mexErrMsgTxt("Error: penalty not supported");

   Matrix<T> duality_gap;
   FISTA::solver<T>(I,convD,alpha0,alpha,param,duality_gap);
   if (nlhs==2) {
      plhs[1]=createMatrix<T>(duality_gap.m(),duality_gap.n());
      T* pr_dualitygap=reinterpret_cast<T*>(mxGetPr(plhs[1]));
      for (int i = 0; i<duality_gap.n()*duality_gap.m(); ++i) pr_dualitygap[i]=duality_gap[i];
   }
   if (param.logName) delete[](param.logName);

}

   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
      if (nrhs != 4)
         mexErrMsgTxt("Bad number of inputs arguments");

      if (nlhs != 1 && nlhs != 2) 
         mexErrMsgTxt("Bad number of output arguments");

      if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
         callFunction<double>(plhs,prhs,nlhs);
      } else {
         callFunction<float>(plhs,prhs,nlhs);
      }
   }




