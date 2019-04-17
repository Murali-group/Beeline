/*!
 * Software SPAMS v2.5 - Copyright 2009-2014 Julien Mairal 
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
 *
 *
 * \file
 *                toolbox lsqsplx
 *
 *                by Yuansi Chen 
 *                yuansi.chen@berkeley.edu
 *
 *                File lsqsplx.h
 * \brief Contains least squares problem with simplex constraint algorithms 
 * It requires the toolbox linalg */

#ifndef LSQSPLX_H
#define LSQSPLX_H

#include <utils.h>
#include "projsplx.h"

/* **************************
 * FISTA: Fast Iterative Shrinkage-Thresholding Algorithm with backtracking
 * **************************/

/// "For" version

template <typename T>
void gpFISTAFor(const Matrix<T>& A, const Vector<T>& b, Vector<T>& xCurr, const T L0 = 1.0, const T eta = 1.0/0.7, const int I = 50, const bool warm = false) {
   const int p = A.n();
   Vector<T> xPrev(p);
   if(warm)
      xPrev.copy(xCurr);
   else
      xPrev.set(1.0/p);
   Vector<T> yPrev(p);
   yPrev.set(1.0/p);
   Vector<T> yCurr(p);

   Vector<T> c(p);
   A.multTrans(b, c, -1.0);

   T L = L0;
   T tPrev = 1.0;
   T tCurr = 1.0;

   Vector<T> gradPrev;
   Vector<T> AyPrev;
   Vector<T> yGrad;
   Vector<T> AyCurr;
   Vector<T> yDelta;
   Vector<T> xDelta;

   for(int i=0; i<I; ++i) {

      A.mult(yPrev, AyPrev);
      A.multTrans(AyPrev, gradPrev);
      gradPrev.add(c);
      T fDelta = 1.0;
      while(true) {

         yGrad.copy(yPrev);
         yGrad.add(gradPrev, -1.0/L);
         projsplx(yGrad, yCurr);
         A.mult(yCurr, AyCurr);

         yDelta.copy(yPrev);
         yDelta.sub(yCurr);
         fDelta = 0.5*AyPrev.nrm2sq() - 0.5*AyCurr.nrm2sq() + c.dot(yDelta);
         if(fDelta + 0.5*L*yDelta.nrm2sq() < gradPrev.dot(yDelta)) {
            L *= eta;
         } else {
            break;
         }
      }

      xCurr.copy(yCurr);
      tCurr = (1.0 + sqrt(1+4*tPrev*tPrev))/2;

      xDelta.copy(xCurr);
      xDelta.sub(xPrev);
      yPrev.copy(xCurr);
      yPrev.add(xDelta, (tPrev-1)/tCurr);
      xPrev.copy(xCurr);
      tPrev = tCurr;
   }
}

/* **************************
 * Active-Set Method with direct inversion, with update(matrix inversion lemma)
 * **************************/

template <typename T>
void activeSet(const Matrix<T>& M, const Vector<T>& b, Vector<T>& xCurr, const T lambda2 = T(1e-5), const T epsilon = T(1e-5), bool warm = false) {
   const int m = M.m();
   const int p = M.n();
   const int L = MIN(m,p)+1;
   Vector<T> c;
   M.multTrans(b, c, -T(1.0));
   T lam2sq =lambda2*lambda2;

   Vector<T> x(p);
   T* const pr_M = M.rawX();
   // constraint matrix
   Vector<T> A(L);
   A.set(T(1.0));
   T* const pr_A = A.rawX();;
   // Non-Active Constraints Set
   Vector<int> NASet(L);
   NASet.set(-1);
   Vector<bool> NAMask(p);
   NAMask.set(false);
   int na;
   Vector<T> xRed(L);
   Vector<T> cRed(L);
   Matrix<T> MRed(m,L);
   T* const pr_MRed = MRed.rawX();
   Matrix<T> GRed(L,L);
   T* const pr_GRed = GRed.rawX();
   Matrix<T> GRedinv(L, L);
   T* const pr_GRedinv = GRedinv.rawX();
   if (!warm) {
      x.setZeros();
      x[0] = T(1.0);
      // Non-Active Constraints Set
      NASet[0] = 0;
      NAMask[0] = true;
      na = 1;
      xRed[0] = x[0];
      cRed[0] = c[0];
      cblas_copy<T>(m, pr_M, 1, pr_MRed, 1);

      // BLAS GRed = MRedT * MRed + lam2sq (na = 1 for now)
      T coeff = cblas_dot<T>(m,pr_MRed,1,pr_MRed,1) + lam2sq;
      GRed(0,0) = coeff;
      GRedinv(0,0)= T(1.0) / GRed(0,0);
   } else {
      // warm start, so use the value of xCurr
      x.copy(xCurr);
      na = 0;
      for(int i = 0; i<p; i++) {
         if(x[i]>T(1e-10)) {
            NASet[na] = i;
            NAMask[i] = true;
            xRed[na] = x[i];
            cRed[na] = c[i];
            cblas_copy<T>(m, pr_M+m*i, 1, pr_MRed+m*na, 1);
            na ++;
         }
      }
      cblas_gemm<T>(CblasColMajor,CblasTrans,CblasNoTrans,na,na,m,T(1.0),pr_MRed,m,pr_MRed,m,T(0.0),
            pr_GRed,L);
      GRed.addDiag(lam2sq);
      GRedinv.copy(GRed);
      sytri<T>(upper,na,pr_GRedinv,L);
      GRedinv.fillSymmetric();
   }

   Vector<T> Gplus(p);
   T* const pr_Gplus = Gplus.rawX();
   Vector<T> Mx;
   M.mult(x, Mx);
   M.multTrans(Mx, Gplus);
   Gplus.add(x, lam2sq);
   Gplus.add(c);

   Vector<T> gRed(L);
   T* const pr_gRed = gRed.rawX();
   Vector<T> MRedxRed(m);

   Vector<T> GinvA(L);
   T* const pr_GinvA = GinvA.rawX();
   Vector<T> Ginvg(L);
   T* const pr_Ginvg = Ginvg.rawX();
   Vector<T> PRed(L);
   T* const pr_PRed = PRed.rawX();
   Vector<T> MRedPRed(m);
   T* const pr_MRedPRed = MRedPRed.rawX();
   Vector<T> UB(L);
   T* const pr_UB = UB.rawX();
   Vector<T> UAiB(L);
   T* const pr_UAiB = UAiB.rawX();
   // main loop active set
   int iter = 0;
   while(iter <= 100*p) {
      ++iter;
      // update of na, NASet, NAMask, xRed, cRed, Gplus, GRedinv, already done
      // now update MRed, gRed, GinvA, Ginvg  (no need to update GRed)
      // MRed
      for(int i = 0; i < na; ++i) {
         // BLAS copy first columns of M to MRed
         cblas_copy<T>(m, pr_M + m*NASet[i], 1, pr_MRed + m*i , 1);
      }
      // gRed
      for(int i = 0; i<na; ++i) {
         gRed[i] = Gplus[NASet[i]];
      }

      // GinvA
      // BLAS GinvA = GRedinv * A (ARed == A)
      cblas_symv<T>(CblasColMajor,CblasUpper,na,T(1.0),pr_GRedinv,L,pr_A,1,T(),pr_GinvA,1);
      // Ginvg
      // BLAS Ginvg = GRedinv * gRed
      cblas_symv<T>(CblasColMajor,CblasUpper,na,T(1.0),pr_GRedinv,L,pr_gRed,1,T(),pr_Ginvg,1);
      T sGinvg = T();
      T sGinvA = T();
      for(int i = 0; i< na; ++i) {
         sGinvg += Ginvg[i];
         sGinvA += GinvA[i];
      }
      T lambdaS = sGinvg / sGinvA;
      // BLAS PRed = GinvA * lambdaS - Ginvg
      cblas_copy<T>(na, pr_GinvA, 1, pr_PRed, 1);
      cblas_scal<T>(na, lambdaS, pr_PRed, 1);
      cblas_axpy<T>(na, T(-1.0), pr_Ginvg, 1, pr_PRed, 1);

      T maxPRed = ABS(PRed[0]);
      for(int i = 0; i< na; ++i) {
         if(ABS(PRed[i])>maxPRed)
            maxPRed = ABS(PRed[i]);
      }
      if(maxPRed < 1e-10) {
         // P = 0, no advance possible
         bool isOpt = true;
         T lamMin = -epsilon;
         int indexMin = -1;
         for(int i = 0; i <p; ++i) {
            if(!NAMask[i] && Gplus[i]-lambdaS < lamMin) {
               isOpt = false;
               lamMin = Gplus[i]-lambdaS;
               indexMin = i;
            }
         }

         if(isOpt) {
            // Got the optimal, STOP!
            xCurr.copy(x);
            return;
         } else {
            // Add one constraint
            NAMask[indexMin] = true;
            NASet[na] = indexMin;
            xRed[na] = x[indexMin];
            cRed[na] = c[indexMin];
            // Gplus inchange

            // update GRedinv
            // BLAS UB = MRed.T * M[:, indexMin]
            cblas_gemv<T>(CblasColMajor,CblasTrans,m,na,T(1.0),pr_MRed,m,pr_M+indexMin*m,1,T(),pr_UB,1);
            // BLAS UC = M[:,indexMin].T* M[:, indexMin]
            T UC = cblas_dot<T>(m,pr_M+indexMin*m,1, pr_M+indexMin*m,1) + lam2sq;
            // BLAS UAiB = GRedinv * UB
            cblas_symv<T>(CblasColMajor,CblasUpper,na,T(1.0),pr_GRedinv,L,pr_UB,1,T(),pr_UAiB,1);
            T USi = T(1.0)/(UC - cblas_dot<T>(na,pr_UB,1,pr_UAiB,1));
            // GRedinv (restricted) += USi * UAiB*UAiB
            //replace cblas_syr<T>(CblasColMajor,CblasUpper,na,USi, pr_UAiB, 1, pr_GRedinv, L);  
            cblas_ger<T>(CblasColMajor,na,na,USi,pr_UAiB,1,pr_UAiB,1,pr_GRedinv,L);
            // copy -UAiB*USi, -UAiB.T*USi, USi to GRedinv
            cblas_copy<T>(na, pr_UAiB, 1, pr_GRedinv+na*L, 1);
            cblas_scal<T>(na, -USi, pr_GRedinv+na*L,1);
            cblas_copy<T>(na, pr_UAiB, 1, pr_GRedinv+na, L);
            cblas_scal<T>(na, -USi, pr_GRedinv+na,L);
            GRedinv(na,na) = USi;

            na += 1;
            assert(na <= L);
         }
      } else {
         // P != 0, can advance
         int indexMin = -1;
         T alphaMin = T(1.0);
         for(int i = 0; i < na; ++i) {
            if(PRed[i] < 0 && -xRed[i]/PRed[i] < alphaMin) {
               indexMin = i;
               alphaMin = -xRed[i]/PRed[i];
            }
         }
         // update x and Gplus
         cblas_scal<T>(na, MIN(T(1.0), alphaMin), pr_PRed, 1);
         for(int i = 0; i< na; ++i) {
            x[NASet[i]] += PRed[i];
            xRed[i] = x[NASet[i]];
            // BLAS Gplus += M.T * M[:, NASet[i]] * cAdv
            //cblas_gemv<T>(CblasColMajor,CblasTrans,m,p,cAdv,pr_M,m, pr_M+NASet[i]*m,1,T(1.0),pr_Gplus,1);
            // BLAS Gplus
            Gplus[NASet[i]] += PRed[i]*lam2sq;
         }
         // Gplus += M.T * MRed * (scaled PRed)
         cblas_gemv<T>(CblasColMajor,CblasNoTrans,m,na,T(1.0),pr_MRed,m,pr_PRed,1,T(),pr_MRedPRed,1);    
         cblas_gemv<T>(CblasColMajor,CblasTrans,m,p,T(1.0),pr_M,m,pr_MRedPRed,1,T(1.0),pr_Gplus,1);


         // delete one constraint or not?
         if(indexMin != -1) {
            // give true 0
            // x[NASet[indexMin]] = T();
            // delete one constraint
            NAMask[NASet[indexMin]] = false;
            // downdate remove this -1;
            na -= 1;
            for(int i = indexMin; i<na; ++i) {
               NASet[i] = NASet[i+1];
               xRed[i] = xRed[i+1];
               cRed[i] = cRed[i+1];
            }
            NASet[na] = -1;
            xRed[na] = T();
            cRed[na] = T();
            // PRed also
            PRed[na] = T();

            // downdate GRedinv
            T UCi = T(1.0)/GRedinv(indexMin, indexMin); 
            // BLAS UB = GRedinv[ALL\indexMin,indexMin]
            cblas_copy<T>(na+1, pr_GRedinv+indexMin*L, 1, pr_UB, 1);
            for(int i = indexMin; i<na; ++i)
               UB[i] = UB[i+1];
            UB[na] = T();
            // get (GRedinv translated)
            // column first
            for(int i = indexMin; i<na; ++i)
               cblas_copy<T>(na+1, pr_GRedinv+(i+1)*L,1,pr_GRedinv+i*L,1);
            // row then
            for(int i = indexMin; i<na; ++i)
               cblas_copy<T>(na+1, pr_GRedinv+i+1, L, pr_GRedinv+i, L);

            // BLAS GRedinv = (GRedinv translated) - UB*UB.T*UCi 
            //replace cblas_syr<T>(CblasColMajor,CblasUpper,na,-UCi, pr_UB, 1, pr_GRedinv, L);  
            cblas_ger<T>(CblasColMajor,na,na,-UCi,pr_UB,1,pr_UB,1,pr_GRedinv,L);
         }
      }
   }
   return;
}

/// Active-Set Method with direct inversion, with update(matrix inversion lemma)
/// Memorize M.T* M + lam2sq = G
template <typename T>
void activeSetS(const Matrix<T>& M, const Vector<T>& b, Vector<T>& xCurr, const Matrix<T>& G, const T lambda2 = 1e-5, const T epsilon = 1e-5, bool warm = false) {
   const int m = M.m();
   const int p = M.n();
   const int L = MIN(m,p)+1;
   T lam2sq =lambda2*lambda2;
   T* const pr_G = G.rawX();
   Vector<T> c;
   M.multTrans(b, c, -T(1.0));

   Vector<T> x(p);
   T* const pr_M = M.rawX();
   // constraint matrix
   Vector<T> A(L);
   A.set(T(1.0));
   T* const pr_A = A.rawX();;
   // Non-Active Constraints Set
   Vector<int> NASet(L);
   NASet.set(-1);
   Vector<bool> NAMask(p);
   NAMask.set(false);
   int na;
   Vector<T> xRed(L);
   Vector<T> cRed(L);
   Matrix<T> MRed(m,L);
   T* const pr_MRed = MRed.rawX();
   Matrix<T> GRed(L,L);
   Matrix<T> GRedinv(L, L);
   T* const pr_GRedinv = GRedinv.rawX();
   Matrix<T> MTMRed(p,L);
   T* const pr_MTMRed = MTMRed.rawX();
   if (!warm) {
      x.setZeros();
      x[0] = T(1.0);
      // Non-Active Constraints Set
      NASet[0] = 0;
      NAMask[0] = true;
      na = 1;
      xRed[0] = x[0];
      cRed[0] = c[0];
      cblas_copy<T>(m, pr_M, 1, pr_MRed, 1);

      // BLAS GRed = MRedT * MRed + lam2sq (na = 1 for now)
      T coeff = cblas_dot<T>(m,pr_MRed,1,pr_MRed,1) + lam2sq;
      GRed(0,0) = coeff;
      GRedinv(0,0)= T(1.0) / GRed(0,0);
      cblas_copy<T>(p, pr_G, 1, pr_MTMRed, 1);
   } else {
      // warm start, so use the value of xCurr
      x.copy(xCurr);
      na = 0;
      for(int i = 0; i<p; i++) {
         if(x[i]>T(1e-10)) {
            NASet[na] = i;
            NAMask[i] = true;
            xRed[na] = x[i];
            cRed[na] = c[i];
            cblas_copy<T>(m, pr_M+m*i, 1, pr_MRed+m*na, 1);
            na ++;
         }
      }
      for(int j = 0; j<na; j++) {
         for(int i = 0; i<na; i++) {
            GRed(j,i) = G(NASet[j],NASet[i]);
         }
      }
      GRedinv.copy(GRed);
      sytri<T>(upper,na,pr_GRedinv,L);
      GRedinv.fillSymmetric();

      for(int j = 0; j<na; j++) {
         cblas_copy<T>(p, pr_G + p*NASet[j], 1, pr_MTMRed + p*j, 1);
      }

   }

   Vector<T> Gplus(p);
   T* const pr_Gplus = Gplus.rawX();
   G.mult(x, Gplus);
   Gplus.add(c);

   Vector<T> gRed(L);
   T* const pr_gRed = gRed.rawX();
   Vector<T> MRedxRed(m);

   Vector<T> GinvA(L);
   T* const pr_GinvA = GinvA.rawX();
   Vector<T> Ginvg(L);
   T* const pr_Ginvg = Ginvg.rawX();
   Vector<T> PRed(L);
   T* const pr_PRed = PRed.rawX();
   Vector<T> UB(L);
   T* const pr_UB = UB.rawX();
   Vector<T> UAiB(L);
   T* const pr_UAiB = UAiB.rawX();
   // main loop active set
   int iter = 0;
   while(iter <= 100*p) {
      ++iter;
      // update of na, NASet, NAMask, xRed, cRed, Gplus, GRedinv, already done
      // now update gRed, GinvA, Ginvg  (no need to update GRed)

      // gRed
      // gRed = Gplus[NASet]
      for(int i = 0; i<na; ++i) {
         gRed[i] = Gplus[NASet[i]];
      }

      // GinvA
      // BLAS GinvA = GRedinv * A (ARed == A)
      cblas_symv<T>(CblasColMajor,CblasUpper,na,T(1.0),pr_GRedinv,L,pr_A,1,T(),pr_GinvA,1);
      // Ginvg
      // BLAS Ginvg = GRedinv * gRed
      cblas_symv<T>(CblasColMajor,CblasUpper,na,T(1.0),pr_GRedinv,L,pr_gRed,1,T(),pr_Ginvg,1);
      T sGinvg = T();
      T sGinvA = T();
      for(int i = 0; i< na; ++i) {
         sGinvg += Ginvg[i];
         sGinvA += GinvA[i];
      }
      T lambdaS = sGinvg / sGinvA;
      // BLAS PRed = GinvA * lambdaS - Ginvg
      cblas_copy<T>(na, pr_GinvA, 1, pr_PRed, 1);
      cblas_scal<T>(na, lambdaS, pr_PRed, 1);
      cblas_axpy<T>(na, T(-1.0), pr_Ginvg, 1, pr_PRed, 1);

      T maxPRed = ABS(PRed[0]);
      for(int i = 0; i< na; ++i) {
         if(ABS(PRed[i])>maxPRed)
            maxPRed = ABS(PRed[i]);
      }
      if(maxPRed < 1e-10) {
         // P = 0, no advance possible
         bool isOpt = true;
         T lamMin = -epsilon;
         int indexMin = -1;
         for(int i = 0; i <p; ++i) {
            if(!NAMask[i] && Gplus[i]-lambdaS < lamMin) {
               isOpt = false;
               lamMin = Gplus[i]-lambdaS;
               indexMin = i;
            }
         }

         if(isOpt) {
            // Got the optimal, STOP!
            xCurr.copy(x);
            return;
         } else {
            // Add one constraint
            NAMask[indexMin] = true;
            NASet[na] = indexMin;
            xRed[na] = x[indexMin];
            cRed[na] = c[indexMin];
            // Gplus inchange
            // update MTMRed
            cblas_copy<T>(p, pr_G+indexMin*p, 1, pr_MTMRed+na*p, 1);

            // update GRedinv
            // BLAS UB = MRed.T * M[:, indexMin]
            // Use G instead here
            for(int i = 0; i< na; ++i) {
               UB[i] = G(NASet[i], indexMin);
            }
            // BLAS UC = M[:,indexMin].T* M[:, indexMin]
            T UC = G(indexMin, indexMin);
            // BLAS UAiB = GRedinv * UB
            cblas_symv<T>(CblasColMajor,CblasUpper,na,T(1.0),pr_GRedinv,L,pr_UB,1,T(),pr_UAiB,1);
            T USi = T(1.0)/(UC - cblas_dot<T>(na,pr_UB,1,pr_UAiB,1));
            // GRedinv (restricted) += USi * UAiB*UAiB
            //replace cblas_syr<T>(CblasColMajor,CblasUpper,na,USi, pr_UAiB, 1, pr_GRedinv, L);  
            cblas_ger<T>(CblasColMajor,na,na,USi,pr_UAiB,1,pr_UAiB,1,pr_GRedinv,L);
            // copy -UAiB*USi, -UAiB.T*USi, USi to GRedinv
            cblas_copy<T>(na, pr_UAiB, 1, pr_GRedinv+na*L, 1);
            cblas_scal<T>(na, -USi, pr_GRedinv+na*L,1);
            cblas_copy<T>(na, pr_UAiB, 1, pr_GRedinv+na, L);
            cblas_scal<T>(na, -USi, pr_GRedinv+na,L);
            GRedinv(na,na) = USi;

            na += 1;
            assert(na <= L);
         }
      } else {
         // P != 0, can advance
         int indexMin = -1;
         T alphaMin = T(1.0);
         for(int i = 0; i < na; ++i) {
            if(PRed[i] < 0 && -xRed[i]/PRed[i] < alphaMin) {
               indexMin = i;
               alphaMin = -xRed[i]/PRed[i];
            }
         }
         // update x and Gplus
         cblas_scal<T>(na, MIN(T(1.0), alphaMin), pr_PRed, 1);
         for(int i = 0; i< na; ++i) {
            x[NASet[i]] += PRed[i];
            xRed[i] = x[NASet[i]];
         }
         // Gplus += MTMRed * (scaled PRed)
         cblas_gemv<T>(CblasColMajor,CblasNoTrans,p,na,T(1.0),pr_MTMRed,p,pr_PRed,1,T(1.0),pr_Gplus,1);


         // delete one constraint or not?
         if(indexMin != -1) {
            // give true 0
            // x[NASet[indexMin]] = T();
            // delete one constraint
            NAMask[NASet[indexMin]] = false;
            // downdate remove this -1;
            na -= 1;
            for(int i = indexMin; i<na; ++i) {
               NASet[i] = NASet[i+1];
               xRed[i] = xRed[i+1];
               cRed[i] = cRed[i+1];
            }
            NASet[na] = -1;
            xRed[na] = T();
            cRed[na] = T();
            // PRed also
            PRed[na] = T();

            // downdate MTMRed
            for(int i = indexMin; i<na; ++i)
               cblas_copy<T>(p, pr_MTMRed+(i+1)*p, 1, pr_MTMRed+i*p, 1);

            // downdate GRedinv
            T UCi = T(1.0)/GRedinv(indexMin, indexMin); 
            // BLAS UB = GRedinv[ALL\indexMin,indexMin]
            cblas_copy<T>(na+1, pr_GRedinv+indexMin*L, 1, pr_UB, 1);
            for(int i = indexMin; i<na; ++i)
               UB[i] = UB[i+1];
            UB[na] = T();
            // get (GRedinv translated)
            // column first
            for(int i = indexMin; i<na; ++i)
               cblas_copy<T>(na+1, pr_GRedinv+(i+1)*L,1,pr_GRedinv+i*L,1);
            // row then
            for(int i = indexMin; i<na; ++i)
               cblas_copy<T>(na+1, pr_GRedinv+i+1, L, pr_GRedinv+i, L);

            // BLAS GRedinv = (GRedinv translated) - UB*UB.T*UCi 
            //replace cblas_syr<T>(CblasColMajor,CblasUpper,na,-UCi, pr_UB, 1, pr_GRedinv, L);  
            cblas_ger<T>(CblasColMajor,na,na,-UCi,pr_UB,1,pr_UB,1,pr_GRedinv,L);
         }
      }
   }
   return;
}

#endif
