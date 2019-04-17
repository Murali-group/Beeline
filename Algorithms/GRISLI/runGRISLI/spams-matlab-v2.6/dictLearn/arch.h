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
 *                toolbox arch
 *
 *                by Yuansi Chen 
 *                yuansi.chen@berkeley.edu
 *
 *                File arch.h
 * \brief Contains archetypal analysis algorithms 
 * It requires the toolbox linalg */

// This file was written by Yuansi Chen
#ifndef ARCH_H
#define ARCH_H

#include <utils.h>
#include "lsqsplx.h"
#include "projsplx.h"

#define NEW_VERSION

/* **************************
 * Alternating Archetypal Analysis 
 * **************************/

/// Alternating Minimization 
/// Each sub-quadratic programming is solved by ActiveSet Method

template <typename T>
void arch(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z,  SpMatrix<T>& A, SpMatrix<T>& B, const int I1 = 3, const int I2 = 20, const T lambda2 = T(10e-5), const T epsilon = T(10e-5),const bool computeZtZ = true);

template <typename T>
void archRobust(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z,  SpMatrix<T>& A, SpMatrix<T>& B, const int I1 = 3, const int I2 = 20, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3),const bool computeZtZ = true);

/// General functions including previous ones. Less parameters and simple use, for Python and Matlab interface

template <typename T>
void archetypalAnalysis(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, SpMatrix<T>& A, SpMatrix<T>& B, const bool robust =false, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50, const int numThreads=-1);

template <typename T>
void archetypalAnalysis(const Matrix<T>& X, Matrix<T>& Z, SpMatrix<T>& A, SpMatrix<T>& B, const bool robust = false, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50, const bool randominit = true, const int numThreads=-1);

template <typename T>
void decompSimplex(const Matrix<T>& X, const Matrix<T>& Z, SpMatrix<T>& alpha, const bool computerZtZ = false, const int numThreads=-1); 

/* **************************
 *  Implementations 
 * **************************/

template <typename T>
void arch(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z,  SpMatrix<T>& A, SpMatrix<T>& B, const int I1, const int I2, const T lambda2, const T epsilon, const bool computeXtX) {
   const int m = X.m();
   const int n = X.n();
   const int p = Z0.n();
   Z.copy(Z0);
   Matrix<T> AlphaT(p,n);
   Matrix<T> BetaT(n,p);
   T RSS = -1.0;
   Vector<T> refColZ;
   Vector<T> copRowAlphaT;
   Vector<T> refColBetaT;
   Matrix<T> matRSS(m,n);
   Vector<T> vBarre(m);
   Vector<T> norms;
   cout.precision(8);

   for(int t=0; t<I1; ++t) {
      // step 1: fix Z to compute Alpha
#pragma omp parallel for
      for(int i=0; i<n; ++i) {
         Vector<T> refColX;
         Vector<T> refColAlphaT;
         X.refCol(i,refColX);
         AlphaT.refCol(i, refColAlphaT);
         gpFISTAFor(Z,refColX, refColAlphaT, T(1.0), T(1.0/0.7), 50, true);
      }
      // step 2: fix Alpha, fix all but one to compute Zi
      Vector<T> refColX;
      for(int l=0; l<p; ++l) {
         AlphaT.copyRow(l, copRowAlphaT);
         T sumAsq =  copRowAlphaT.nrm2sq();
         Z.refCol(l, refColZ);
         // matRSS = X- Z*AlphaT
         matRSS.copy(X);
         Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
         if(sumAsq < T(10e-8)) {
            // singular
            matRSS.norm_2_cols(norms);
            int k = norms.max();
            X.refCol(k, refColX);
            refColZ.copy(refColX);
         } else {
            matRSS.rank1Update(refColZ, copRowAlphaT);
            matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
            // least square to get Beta
            BetaT.refCol(l, refColBetaT);
            gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 50, true);
            X.mult(refColBetaT, refColZ);
         }
      }

      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      RSS = matRSS.normFsq();
      cout << "RSS FISTA = " << RSS << endl;
      flush(cout);
   }  

   for(int t=0; t<I2; ++t) {
      Matrix<T> G;
      if (computeXtX) {
         Z.XtX(G);
         G.addDiag(lambda2*lambda2);
      }
      // step 1: fix Z to compute Alpha
#pragma omp parallel for
      for(int i=0; i<n; ++i) {
         Vector<T> refColX;
         Vector<T> refColAlphaT;
         X.refCol(i,refColX);
         AlphaT.refCol(i, refColAlphaT);
         if (computeXtX) {
            activeSetS<T>(Z,refColX, refColAlphaT, G, lambda2, epsilon);
         } else {
            activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon);
         }
      }
      // step 2: fix Alpha, fix all but one to compute Zi
#ifdef NEW_VERSION
      // new version
      Vector<T> refColX;
      Vector<T> tmp;
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      for(int l=0; l<p; ++l) {
         AlphaT.copyRow(l, copRowAlphaT);
         T sumAsq =  copRowAlphaT.nrm2sq();
         Z.refCol(l, refColZ);
         // matRSS = X- Z*AlphaT
         if(sumAsq < T(10e-8)) {
            // singular
            matRSS.norm_2_cols(norms);
            int k = norms.max();
            X.refCol(k, refColX);
            refColZ.copy(refColX);
         } else {
            //matRSS.rank1Update(refColZ, copRowAlphaT);
            matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
            vBarre.add(refColZ);
            tmp.copy(refColZ);
            // least square to get Beta
            BetaT.refCol(l, refColBetaT);
            activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
            X.mult(refColBetaT, refColZ);
            tmp.sub(refColZ);
            matRSS.rank1Update(tmp, copRowAlphaT);
         }
      }
#else
      // end new version

      Vector<T> refColX;
      for(int l=0; l<p; ++l) {
         AlphaT.copyRow(l, copRowAlphaT);
         T sumAsq =  copRowAlphaT.nrm2sq();
         Z.refCol(l, refColZ);
         // matRSS = X- Z*AlphaT
         matRSS.copy(X);
         Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
         if(sumAsq < T(10e-8)) {
            // singular
            matRSS.norm_2_cols(norms);
            int k = norms.max();
            X.refCol(k, refColX);
            refColZ.copy(refColX);
         } else {
            matRSS.rank1Update(refColZ, copRowAlphaT);
            matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
            // least square to get Beta
            BetaT.refCol(l, refColBetaT);
            activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
            X.mult(refColBetaT, refColZ);
         }
      }

      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
#endif
      RSS = matRSS.normFsq();
      cout << "RSS AS = " << RSS << endl;
      flush(cout);
   }
   AlphaT.toSparse(A);
   BetaT.toSparse(B);
}

template <typename T>
void archRobust(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z,  SpMatrix<T>& A, SpMatrix<T>& B, const int I1, const int I2, const T lambda2, const T epsilon, const T epsilon2, const bool computeXtX) {
   const int m = X.m();
   const int n = X.n();
   const int p = Z0.n();
   Z.copy(Z0);
   Matrix<T> AlphaT(p,n);
   Matrix<T> BetaT(n,p);

   T RSN = -1.0;
   Vector<T> refColZ;
   Vector<T> copRowAlphaT;
   Vector<T> refColBetaT;
   Matrix<T> matRSS(m,n);
   Vector<T> vBarre(m);
   Vector<T> norms;
   cout.precision(8);

   for(int t=0; t<I1; ++t) {
      // step 1: fix Z to compute Alpha
#pragma omp parallel for
      for(int i=0; i<n; ++i) {
         Vector<T> refColX;
         Vector<T> refColAlphaT;
         X.refCol(i,refColX);
         AlphaT.refCol(i, refColAlphaT);
         gpFISTAFor(Z, refColX, refColAlphaT, T(1.0), T(1.0/0.7), 10, true);
      }
      // update scale factors
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      matRSS.norm_2_cols(norms);
      norms.thrsmax(epsilon2);
      norms.Sqrt();
      Vector<T> refColX;
      // step 2: fix Alpha, fix all but one to compute Zi
      for(int l=0; l<p; ++l) {
         Z.refCol(l, refColZ);

         AlphaT.copyRow(l, copRowAlphaT);
         copRowAlphaT.div(norms);
         T sumAsq =  copRowAlphaT.nrm2sq();

         matRSS.copy(X);
         Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

         if(sumAsq < T(10e-8)) {
            // singular
            matRSS.norm_2_cols(norms);
            int k = norms.max();
            X.refCol(k, refColX);
            refColZ.copy(refColX);
         } else {
            // absorbe the weights by rowAlphaT
            copRowAlphaT.div(norms);
            matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
            vBarre.add(refColZ);
            // least square to get Beta
            BetaT.refCol(l, refColBetaT);
            gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 10, true); 
            X.mult(refColBetaT, refColZ);
         }
      }

      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      matRSS.norm_2_cols(norms);
      for (int i=0; i<norms.n(); ++i)
         if (norms[i] <= epsilon2)
            norms[i]=norms[i]*norms[i]/(2*epsilon2) + epsilon2/2;
      RSN = norms.sum();
      cout << "RSN FISTA= " << RSN << endl;
      flush(cout);
   }

   for(int t=0; t<I2; ++t) {
      Matrix<T> G;
      if (computeXtX) {
         Z.XtX(G);
         G.addDiag(lambda2*lambda2);
      }
      // step 1: fix Z to compute Alpha
#pragma omp parallel for
      for(int i=0; i<n; ++i) {
         Vector<T> refColX;
         Vector<T> refColAlphaT;
         X.refCol(i,refColX);
         AlphaT.refCol(i, refColAlphaT);
         if (computeXtX) {
            activeSetS<T>(Z,refColX, refColAlphaT, G, lambda2, epsilon);
         } else {
            activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon);
         }
      }
      // update scale factors
#ifndef NEW_VERSION
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      matRSS.norm_2_cols(norms);
      norms.thrsmax(epsilon2);
      norms.Sqrt();
      // step 2: fix Alpha, fix all but one to compute Zi
      Vector<T> refColX;
      for(int l=0; l<p; ++l) {
         Z.refCol(l, refColZ);

         AlphaT.copyRow(l, copRowAlphaT);
         copRowAlphaT.div(norms);
         T sumAsq =  copRowAlphaT.nrm2sq();

         matRSS.copy(X);
         Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

         if(sumAsq < T(10e-8)) {
            // singular
            matRSS.norm_2_cols(norms);
            int k = norms.max();
            X.refCol(k, refColX);
            refColZ.copy(refColX);
         } else {
            // absorbe the weights by rowAlphaT
            copRowAlphaT.div(norms);
            matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
            vBarre.add(refColZ);
            // least square to get Beta
            BetaT.refCol(l, refColBetaT);
            activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
            X.mult(refColBetaT, refColZ);
         }
      }

      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
#else
      /// new version
      Vector<T> refColX;
      Vector<T> tmp;
      Vector<T> tmp2;
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      matRSS.norm_2_cols(norms);
      norms.thrsmax(epsilon2);
      norms.Sqrt();
      for(int l=0; l<p; ++l) {
         Z.refCol(l, refColZ);

         AlphaT.copyRow(l, copRowAlphaT);
         tmp2.copy(copRowAlphaT);
         copRowAlphaT.div(norms);
         T sumAsq =  copRowAlphaT.nrm2sq();

         if(sumAsq < T(10e-8)) {
            // singular
            matRSS.norm_2_cols(tmp);
            int k = tmp.max();
            X.refCol(k, refColX);
            refColZ.copy(refColX);
         } else {
            // absorbe the weights by rowAlphaT
            copRowAlphaT.div(norms);
            matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
            vBarre.add(refColZ);
            tmp.copy(refColZ);
            // least square to get Beta
            BetaT.refCol(l, refColBetaT);
            activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
            X.mult(refColBetaT, refColZ);
            tmp.sub(refColZ);
            matRSS.rank1Update(tmp,tmp2);
         }
      }
#endif
      /// end new version
      
      matRSS.norm_2_cols(norms);
      for (int i=0; i<norms.n(); ++i)
         if (norms[i] <= epsilon2)
            norms[i]=norms[i]*norms[i]/(2*epsilon2) + epsilon2/2;
      RSN = norms.sum();
      cout << "RSN AS= " << RSN << endl;
      flush(cout);
   }
   AlphaT.toSparse(A);
   BetaT.toSparse(B);
}

template <typename T>
void archetypalAnalysis(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, SpMatrix<T>& A, SpMatrix<T>& B, const bool robust, const T epsilon2, const bool computeXtX, const int stepsFISTA, const int stepsAS, const int numThreads) {
   init_omp(numThreads);
   const T epsilon = 1e-5;
   const T lambda2 = 1e-5;
   if (!robust) {
      arch(X, Z0, Z, A, B, stepsFISTA, stepsAS, epsilon,lambda2,computeXtX);
   } else {
      archRobust(X, Z0, Z, A, B, stepsFISTA, stepsAS, epsilon,lambda2,epsilon2,computeXtX);
   }
}

template <typename T>
void archetypalAnalysis(const Matrix<T>& X, Matrix<T>& Z, SpMatrix<T>& A, SpMatrix<T>& B, const bool robust, const T epsilon2, const bool computeXtX, const int stepsFISTA, const int stepsAS, const bool randominit, const int numThreads) {

   const int m = X.m();
   const int n = X.n();
   const int p = Z.n();
   Matrix<T> Z0(m,p);
   Vector<T> refColZ0;
   Vector<T> refColX;
   if(!randominit) {
      for(int i=0; i<p; i++) {
         X.refCol(i%n, refColX);
         Z0.refCol(i%n, refColZ0);
         refColZ0.copy(refColX);
      }
   } else {
      srandom(0);
      for(int i=0; i<p; i++) {
         int k = random() % n;
         X.refCol(k, refColX);
         Z0.refCol(i, refColZ0);
         refColZ0.copy(refColX);
      }
   }
   archetypalAnalysis(X, Z0, Z, A, B, robust, epsilon2, computeXtX, stepsFISTA, stepsAS,numThreads);
}

template <typename T>
void decompSimplex(const Matrix<T>& X, const Matrix<T>& Z, SpMatrix<T>& alpha, const bool computeZtZ, const int numThreads) {
   init_omp(numThreads);
   const int n = X.n();
   const int p = Z.n();
   Matrix<T> AlphaT(p,n);
   int i;
   if(computeZtZ) {
      Matrix<T> G;
      Z.XtX(G);
      T lambda2 = 1e-5;
      G.addDiag(lambda2*lambda2);
#pragma omp parallel for private(i)
      for(i=0; i<n; ++i) {
         Vector<T> refColX;
         Vector<T> refColAlphaT;
         X.refCol(i,refColX);
         AlphaT.refCol(i, refColAlphaT);
         activeSetS(Z,refColX, refColAlphaT, G);
      }
      AlphaT.toSparse(alpha);
   } else {
#pragma omp parallel for private(i)
      for(i=0; i<n; ++i) {
         Vector<T> refColX;
         Vector<T> refColAlphaT;
         X.refCol(i,refColX);
         AlphaT.refCol(i, refColAlphaT);
         activeSet(Z,refColX, refColAlphaT);
      }
      AlphaT.toSparse(alpha);
   }
}

#endif
