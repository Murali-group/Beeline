#ifndef SVM_H
#define SVM_H  

#include <linalg.h>

template <typename T>
T normsq(const Vector<T>& x, const Vector<T>& y) {
   return x.nrm2sq()+y.nrm2sq()-2*y.dot(x);
}


template <typename T>
void miso_svm_aux(const Vector<T>& y, const Matrix<T>& X, Vector<T>& w, const T R, const T lambda, const T eps, const int max_iter) {
   const int n = y.n();
   w.setZeros();
   const T L = R+lambda;
   const T deltaT = n*MIN(T(1.0)/n,lambda/(2*L));
   Vector<T> xi;
   Vector<T> alpha(n);
   alpha.setZeros();
   Vector<T> C(n);
   C.setZeros();
   Vector<T> tmp;
   T dualold=0;
   T dual=0;
   for (int ii = 0; ii<max_iter; ++ii) {
      if (ii > 0 && (ii % (10*n)) == 0) {
         X.mult(alpha,w,T(1.0)/n); // to improve numerical stability
         X.multTrans(w,tmp);
         T primal=0;
         for (int kk=0; kk<n; ++kk) {
            const T los=MAX(0,1-y[kk]*tmp[kk]);
            primal += los*los;
         }
         primal *= T(0.5)/n;
         T reg=0.5*lambda*w.nrm2sq();
         primal += reg;
         dual=C.mean() - reg;
         if (dual <= dualold || (primal - dual) < eps) {
#pragma omp critical
            {
               cout << "Solver has finished after " << ii << " iterations, primal: " << primal << ", dual: " << dual << ", gap: " << (primal-dual) << endl;
            }
            break;
         }
         dualold=dual;
      }
      const int ind = random() % n;
      const T yi=y[ind];
      X.refCol(ind,xi);
      const T beta = yi*xi.dot(w);
      const T gamma=MAX(T(1.0)-beta,0);
      T newalpha;
      C[ind]=(T(1.0)-deltaT)*C[ind]+deltaT*(T(0.5)*gamma*gamma+beta*gamma);
      newalpha=(T(1.0)-deltaT)*(alpha[ind])+deltaT*yi*gamma/lambda;
      w.add(xi,(newalpha-alpha[ind])/n);
      alpha[ind]=newalpha;
   }
};

template <typename T>
void miso_svm_onevsrest(const Vector<T>& yAll, const Matrix<T>& X, Matrix<T>& W, const T lambda, const T eps, const int max_iter, const bool accelerated = false) {
   const int n = yAll.n();
   const int p = X.m();
   const int nclasses=yAll.maxval()+1;
   W.resize(p,nclasses);
   Vector<T> normX;
   X.norm_2sq_cols(normX);
   const T R = normX.mean();
   cout << "Value of R: " << R << endl;

   cout << "Problem size: p x n: " << p << " " << n << endl;
   cout << "*********************" << endl;
   cout << "Processes Lambda " << lambda << endl;
   cout << "Eps " << eps << endl;
   int jj;
#pragma omp parallel for private(jj)
   for (jj = 0; jj<nclasses; ++jj) {
      Vector<T> w;
      W.refCol(jj,w);
      Vector<T> y(n);
      for (int ii = 0; ii<n; ++ii) 
         y[ii]= abs<T>((yAll[ii] - T(jj))) < T(0.1) ? T(1.0) : -T(1.0);
      if (accelerated && T(2.0)*R/n > lambda) {
         accelerated_miso_svm_aux(y,X,w,R,lambda,eps,max_iter);
      } else {
         miso_svm_aux(y,X,w,R,lambda,eps,max_iter);
      }
   }
}

template <typename T>
void miso_svm(const Vector<T>& y, const Matrix<T>& X, Matrix<T>& W, const Vector<T>& tablambda, const T eps, const int max_iter) {
   const int n = y.n();
   const int p = X.m();
   const int nlambda=tablambda.n();
   W.resize(p,nlambda);
   W.setZeros();
   Vector<T> normX;
   X.norm_2sq_cols(normX);
   const T R = normX.fmax();

   cout << "Problem size: p x n: " << p << " " << n << endl;
   for (int jj = 0; jj<nlambda; ++jj) {
      const T lambda=tablambda[jj];
      cout << "*********************" << endl;
      cout << "Processes Lambda " << lambda << endl;
      Vector<T> w;
      W.refCol(jj,w);
      miso_svm_aux(y,X,w,R,lambda,eps,max_iter);
   }
}

template <typename T>
void accelerated_miso_svm_aux(const Vector<T>& y, const Matrix<T>& X, Vector<T>& w, const T R, const T lambda, const T eps, const int max_iter) {
   const int n = y.n();
   const int p = X.m();
   w.setZeros();
   Vector<T> alpha(n);
   alpha.setZeros();
   Vector<T> C(n);
   C.setZeros();
   Vector<T> z(p);
   z.setZeros();
   Vector<T> zold(p);
   zold.setZeros();
   Vector<T> wold(p);
   wold.setZeros();
   Vector<T> xtw(n);
   xtw.setZeros();
   const T kappa = (T(2.0)*R/n-lambda); 
   const T q = lambda/(lambda+kappa);
   const T qp = T(0.9)*sqrt(q);
   const T alphak = sqrt(q);
   const T betak=(T(1.0)-alphak)/(T(1.0)+alphak);
   T epsk=T(1.0);
   T gapk=T(1.0);
   T gap=T(1.0);
   int total_iters=0;
   int counter = 1;
   T gapold=T(1.0);
   for (int ii=0; ii<max_iter; ++ii) {
      epsk *= (T(1.0)-qp);
      // check if continue or not
      wold.copy(w);
      if ((total_iters / (10*n)) >= counter) {
         ++counter;
         w.copy(z);
         w.scal(kappa/(kappa+lambda));
         X.mult(alpha,w,lambda/(n*(kappa+lambda)),T(1.0));
      } else {
         w.add(z,kappa/(kappa+lambda));
         w.add(zold,-kappa/(kappa+lambda));
      }
      const T diffNorm = normsq(z,zold);
      gapk=(n*(gapk + T(0.5)*(kappa*kappa/(lambda+kappa))*diffNorm));
      T loss;
      int num_iters;
      accelerated_miso_svm_aux2(y, X, w, alpha, C, loss, gapk, num_iters, z, kappa, R, lambda, epsk);
      total_iters += num_iters;
      const T primal = loss+T(0.5)*lambda*w.nrm2sq();
      Vector<T> ws;
      ws.copy(w);
      ws.scal((kappa+lambda)/lambda);
      ws.add(z,-kappa/lambda);
      const T dual=C.mean() - T(0.5)*lambda*ws.nrm2sq();
      gap=primal-dual;
      if ((ii > 30 && gap >= gapold) || gap <= eps || total_iters >= max_iter) {
#pragma omp critical
         {
            cout << "Iteration " << total_iters << ", inner it: " << ii << ", loss: " << loss << ", primal: " << primal << ", dual: " << dual << ", gap: " << (primal-dual) << endl;
         }
         break;
      }
      gapold=gap;
      zold.copy(z);
      z.copy(w);
      z.scal(T(1.0)+betak);
      z.add(wold,-betak);
   }
};


// need to restart !
template <typename T>
void accelerated_miso_svm_aux2(const Vector<T>& y, const Matrix<T>& X, Vector<T>& w, Vector<T>& alpha, Vector<T>& C, T& loss,T& gap, int& num_iters,  const Vector<T>& z, const T kappa, const T R, const T lambda, const T eps) {
   const int n = y.n();
//   const int p = X.m();
   const long long max_iter = static_cast<long long>(floor(log(double(eps)/double(gap))/log(double(1.0)-double(1.0)/n)));
   Vector<T> tmp;
   Vector<T> xi;
   //const T deltaT = n*MIN(T(1.0)/n,(lambda+kappa)/(2*R));
   for (int ii = 0; ii<max_iter; ++ii) {
      if (ii > 0 && (ii % (n)) == 0) {
         loss=0;
         X.multTrans(w,tmp);
         for (int kk=0; kk<n; ++kk) {
            const T los=MAX(0,1-y[kk]*tmp[kk]);
            loss += los*los;
         }
         loss *= T(0.5)/n;
         const T reg=T(0.5)*(lambda+kappa)*w.nrm2sq();
         const T primal = loss+ reg  - kappa*w.dot(z);
         const T dual=C.mean() - reg;
         if ((primal - dual) < eps) {
   //         cout << "   Inner iteration " << ii << ", loss: " << loss << ", primal: " << primal << ", dual: " << dual << ", gap: " << (primal-dual) << "eps: " << eps << endl;
            gap=primal-dual;
            num_iters=ii;
            break;
         }
      }
      const int ind = random() % n;
      const T yi=y[ind];
      X.refCol(ind,xi);
      const T beta = yi*xi.dot(w);
      const T gamma=MAX(T(1.0)-beta,0);
      T newalpha;
      //if (deltaT > 0.999999999999999999) {
         C[ind]=T(0.5)*gamma*gamma+beta*gamma;
         newalpha=yi*gamma/lambda;
         if (newalpha != alpha[ind])
            w.add(xi,lambda*(newalpha-alpha[ind])/(n*(lambda+kappa)));
      //} else {
      //   C[ind]=(T(1.0)-deltaT)*C[ind]+deltaT*(T(0.5)*gamma*gamma+beta*gamma);
      //   newalpha=(T(1.0)-deltaT)*(alpha[ind])+deltaT*yi*gamma/lambda;
      //   w.add(xi,lambda*(newalpha-alpha[ind])/(n*(lambda+kappa)));
     //` }
      alpha[ind]=newalpha;
   }
}


#endif
