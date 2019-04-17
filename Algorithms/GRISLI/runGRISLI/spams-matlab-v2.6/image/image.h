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


#ifndef IMAGE_H
#define IMAGE_H

#include <linalg.h>

void buildGrid(Vector<INTM>& px, Vector<INTM>& py, const INTM nx,
      const INTM ny, const INTM n, const INTM st);

/// class Image
template<typename T> class Image : public Data<T> {
   public:
      /// Constructor_ind
      Image(INTM w, INTM h, INTM numChannels = 1);
      /// Constructor with existing data 
      Image(T* X, INTM w, INTM h, INTM numChannels = 1);
      /// Empty constructor
      Image();

      /// Destructor
      ~Image();

      /// clear the image
      inline void clear();

      /// Accessor
      inline INTM numChannels() const  { return _V;};
      /// Accessor
      inline INTM width() const  { return _w;};
      /// Accessor
      inline INTM height() const  { return _h;};
      /// Low-level accessor (DANGEROUS)
      inline T* rawX() const { return _X;};
      /// operator (i,j)
      inline T operator[](const INTM i) const { return _X[i];};
      /// operator (i,j)
      inline T& operator()(const INTM i) const { return _X[i];};
      /// operator (i,j)
      inline T& operator()(const INTM i, const INTM j) const { return _X[j*_h+i];};
      /// operator (i,j)
      inline T& operator()(const INTM i, const INTM j, const INTM k) const { return _X[_h*_w*k+j*_h+i];};
      /// accessor
      inline INTM n() const { return _n;};
      /// accessor
      inline INTM m() const { return _m;};
      /// accessor
      inline INTM V() const { return _V;};

      /// Extract the patches of an image
      void extractPatches(Matrix<T>& X, const INTM n, 
            const INTM st = 1); 
      void combinePatches(const Matrix<T>& X, const T lambda = 0, const INTM st = 1, const bool normalize=true); 
      void getPatch(Vector<T>& patch, const INTM i, const INTM n) const;
      void getPatch(Vector<T>& patch, const INTM i, const INTM j, const INTM n) const;
      void getData(Vector<T>& data, const INTM i) const;
      void getGroup(Matrix<T>& data, const vector_groups& groups,
            const INTM i) const;
      inline void scal(const T lambda);
      inline void copy(const Image<T>& image);
      inline void setZeros();
      inline void resize(const INTM w, const INTM h, const INTM V = 1);

   private:
      /// Forbid lazy copies
      explicit Image<T>(const Image<T>& image);
      /// Forbid lazy copies
      Image<T>& operator=(const Image<T>& image);

      /// is the data allocation external or not
      bool _externAlloc;
      /// poINTMer to the data
      T* _X;
      /// number of rows
      INTM _h;
      /// number of columns
      INTM _w;
      /// number of channels
      INTM _V;

      /// data variable
      INTM _n;
      /// data variable
      INTM _m;
      /// size of patch
      INTM _sizePatch;
};

/// Empty constructor
template <typename T> Image<T>::Image() :
 _externAlloc(false), _X(NULL), _h(0), _w(0), _V(1), _n(0),
_m(0),_sizePatch(0) { };

/// Destructor
template <typename T> Image<T>::~Image() {
   clear();
}

/// Constructor with existing data X of an image
template <typename T> Image<T>::Image(INTM w, INTM h, INTM numChannels) :
   _externAlloc(false), _h(h) ,_w(w), _V(numChannels),_n(0),_m(0),_sizePatch(0) { 
      _X=new T[_w*_h*_V];
   };

/// Constructor with existing data X of an image
template <typename T> Image<T>::Image(T* X, INTM w, INTM h, INTM numChannels) :
   _externAlloc(true), _X(X), _h(h), _w(w), _V(numChannels),_n(0),_m(0),_sizePatch(0) {  };

/// Clear the image
template <typename T> inline void Image<T>::clear() {
   if (!_externAlloc) delete[](_X);
   _w=0;
   _h=0;
   _V=1;
   _X=NULL;
   _externAlloc=true;
};

/// Extract all the patches of an image
template <typename T>
void Image<T>::extractPatches(Matrix<T>& Xm, 
      const INTM n, const INTM st) {
   T* I = _X;
   const INTM ny = _w;
   const INTM nx = _h;
   const INTM V = _V;
   const INTM sizePatch=n*n*V;
   const INTM sizeChannel=n*n;
   Vector<INTM> px, py;
   buildGrid(px,py,nx,ny,n,st);
   const INTM numPatchesX=px.n();
   const INTM numPatchesY=py.n();
   const INTM numPatches=px.n()*py.n();
   Xm.resize(n*n*V,numPatches);
   T* X = Xm.rawX();
   INTM* x = px.rawX();
   INTM* y = py.rawX();
   const INTM nxy = ny*nx;
   for (INTM i = 0; i<V; ++i) {
      INTM count=0;
      for (INTM jj = 0; jj<numPatchesY; ++jj) {
         const INTM j=y[jj];
         for (INTM kk = 0; kk<numPatchesX; ++kk) {
            const INTM k = x[kk];
            for (INTM l = 0; l<n; ++l) { 
               memcpy(X+count*sizePatch+i*sizeChannel+n*l,I+i*nxy+(j+l)*nx+k,n*sizeof(T));
            } 
            ++count;
         }
      }
   }
};

/// Extract all the patches of an image
template <typename T>
void Image<T>::combinePatches(const Matrix<T>& Xm, const T lambda,
      const INTM st, const bool normalize) {
   const INTM n = static_cast<INTM>(sqrt(static_cast<T>(Xm.m()/_V)));
   if (Xm.m() != n*n*_V)
      return;
   T* X = Xm.rawX();
   //this->setZeros();
   T* I = _X;
   const INTM ny = _w;
   const INTM nx = _h;
   const INTM V = _V;
   const INTM sizePatch=n*n*V;
   const INTM sizeChannel=n*n;
   Vector<INTM> px, py;
   buildGrid(px,py,nx,ny,n,st);
   const INTM numPatchesX=px.n();
   const INTM numPatchesY=py.n();
   const INTM numPatches=px.n()*py.n();
   INTM* x = px.rawX();
   INTM* y = py.rawX();

   this->scal(lambda);

   if (sizePatch != Xm.m() || numPatches != Xm.n()) {
      return;
   }

   const INTM nxy = ny*nx;
   for (INTM i = 0; i<V; ++i) {
      INTM count=0;
      for (INTM jj = 0; jj<numPatchesY; ++jj) {
         const INTM j=y[jj];
         for (INTM kk = 0; kk<numPatchesX; ++kk) {
            const INTM k = x[kk];
            for (INTM l = 0; l<n; ++l) { 
               cblas_axpy<T>(n,T(1.0),X+count*sizePatch+i*sizeChannel+n*l,1,
                     I+i*nxy+(j+l)*nx+k,1); 
            } 
            ++count;
         }
      }
   }

   if (normalize) {
      T* den = new T[_w*_h];
      for (INTM i = 0; i<_w*_h; ++i)
         den[i]=lambda;
      INTM count=0;
      for (INTM jj = 0; jj<numPatchesY; ++jj) {
         const INTM j=y[jj];
         for (INTM kk = 0; kk<numPatchesX; ++kk) {
            const INTM k = x[kk];
            for (INTM l = 0; l<n; ++l) 
               for (INTM m = 0; m<n; ++m) den[(j+l)*nx+k+m]++; 
            ++count;
         }
      }

      for (INTM i = 0; i<V; ++i)
         vDiv<T>(_w*_h,I+i*nxy,den,I+i*nxy);
      delete[](den);
   }
};

template <typename T>
inline void Image<T>::copy(const Image<T>& image) {
   const INTM w=image.width();
   const INTM h=image.height();
   const INTM V=image.numChannels();
   this->resize(w,h,V);
   memcpy(_X,image.rawX(),_w*_h*_V*sizeof(T));
};

/// resize the image
template <typename T>
inline void Image<T>::resize(const INTM w, const INTM h, const INTM V) {
   INTM actSize=_w*_h*_V;
   if (actSize == w*h*V) {
      _w = w;
      _h = h;
      _V = V;
      return;
   }
   clear();
   _w = w;
   _h = h;
   _V = V;
   _externAlloc=false;
   if (_w*_h*_V == 0) {
      _X = NULL;
      return;
   }
   _X=new T[_w*_h*_V];
};

template <typename T>
inline void Image<T>::setZeros() {
   Vector<T> X(_X,_h*_w*_V);
   X.setZeros();
};

void buildGrid(Vector<INTM>& px, Vector<INTM>& py, const INTM nx,
      const INTM ny, const INTM n, const INTM st) {
   const INTM numPatchesX = st == 1 ? nx-n+1 : (nx-n+1)/st + 1;
   const INTM numPatchesY = st == 1 ? ny-n+1 : (ny-n+1)/st + 1;
   px.resize(numPatchesX);
   py.resize(numPatchesY);
   INTM * prx=px.rawX();
   INTM * pry=py.rawX();
   INTM count=0;
   for (INTM j = 0; j<st*(numPatchesY-1); j += st) 
      pry[count++]=j;
   pry[count]=ny-n;
   count=0;
   for (INTM k = 0; k<st*(numPatchesX-1); k += st)  {
      prx[count++]=k;
   }
   prx[count]=nx-n;
};


template <typename T>
inline void Image<T>::getPatch(Vector<T>& patch, const INTM i, const INTM n) const  {
   const INTM jj = i/(_h-n+1);
   const INTM ii = i-jj*(_h-n+1);
   this->getPatch(patch,ii,jj,n);
};

template <typename T>
inline void Image<T>::getData(Vector<T>& data, const INTM i) const {
   this->getPatch(data,i,_sizePatch);
};

template <typename T>
inline void Image<T>::getGroup(Matrix<T>& data, const vector_groups& groups,
      const INTM i) const {
   const group& gr = groups[i];
   const INTM N = gr.size();
   data.resize(_sizePatch*_sizePatch*_V,N);
   Vector<T> d;
   INTM count=0;
   for (group::const_iterator it = gr.begin(); it != gr.end(); ++it) {
      data.refCol(count,d);
      this->getPatch(d,*it,_sizePatch);
      ++count;
   }
};

template <typename T>
inline void Image<T>::getPatch(Vector<T>& patch, const INTM i, const INTM j,
      const INTM n) const {
   T* patchX = patch.rawX();
   for (INTM ii = 0; ii<_V; ++ii) {
      for (INTM jj = 0;jj<n; ++jj) {
         //cblas_copy<T>(n,_X+ii*_h*_w+(jj+j)*_h+i,1,patchX+ii*n*n+jj*n,1);
         memcpy(patchX+ii*n*n+jj*n,_X+ii*_h*_w+(jj+j)*_h+i,n*sizeof(T));
      }
   }
};

template <typename T>
inline void Image<T>::scal(const T lambda) {
   Vector<T> tmp(_X,_h*_w*_V);
   tmp.scal(lambda);
}

template <typename T> class ConvolveDictionary : public AbstractMatrixB<T> {
   public:
      inline INTM n() const { return _m; };
      inline INTM m() const { return _p; };

      ConvolveDictionary(const Matrix<T>& D, const INTM w, const INTM h, const INTM nim = 1);

      inline void print(const string& name) const;

      virtual T dot(const Matrix<T>& x) const {
      FLAG(45)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void copyTo(Matrix<T>& copy) const {
      FLAG(46)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void copyRow(const INTM i, Vector<T>& x) const {
      FLAG(47)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void XtX(Matrix<T>& XtX) const {
      FLAG(48)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void multSwitch(const Matrix<T>& B, Matrix<T>& C, 
            const bool transA = false, const bool transB = false,
            const T a = 1.0, const T b = 0.0) const {
      FLAG(49)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void mult(const Matrix<T>& B, Matrix<T>& C, 
            const bool transA = false, const bool transB = false,
            const T a = 1.0, const T b = 0.0) const {
      FLAG(50)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void mult(const SpMatrix<T>& B, Matrix<T>& C, 
            const bool transA = false, const bool transB = false,
            const T a = 1.0, const T b = 0.0) const {
      FLAG(51)
         std::cerr << "Not Implemented" << std::endl;
         exit(1);
      }
      virtual void multTrans(const Vector<T>& x, Vector<T>& b,
            const T alpha = 1.0, const T beta = 0.0) const;

      /// perform b = alpha*A*x + beta*b, when x is sparse
      virtual void mult(const SpVector<T>& x, Vector<T>& b, 
            const T alpha = 1.0, const T beta = 0.0) const;

      virtual void mult(const Vector<T>& x, Vector<T>& b, 
            const T alpha = 1.0, const T beta = 0.0) const;

   private:
      INTM _w; // width of images
      INTM _h; // height of images
      INTM _nim;
      INTM _sizeEdges;  // size of patches
      INTM _m;  // number of pixels, virtual size of dictionary elements
      INTM _p;  // virtual number of dictionary elements
      INTM _sizeMapX;
      INTM _sizeMapY;
      Matrix<T> _D;
};

template <typename T>
ConvolveDictionary<T>::ConvolveDictionary(const Matrix<T>& D, const INTM w, const INTM h, const INTM nim) {
   _w=w;
   _h=h;
   _nim=nim;
   _D.copy(D);
   _sizeEdges=static_cast<INTM>(sqrt(static_cast<T>(D.m()/nim)));
   _m=_h*_w*_nim;
   _sizeMapX=_h-_sizeEdges+1;
   _sizeMapY=_w-_sizeEdges+1;
   _p=D.n()*_sizeMapX*_sizeMapY;
};

template <typename T>
void ConvolveDictionary<T>::print(const string& name) const {
   std::cerr << "Convolutional Dictionary" << std::endl;
   std::cerr << name << std::endl;
   std::cerr << _m << " x " << _p << std::endl;
   std::cerr << "Image: " << _w << " x " << _h << " x " << _nim << std::endl;
   std::cerr << "Patches: " << _sizeEdges << " x " << _sizeEdges << std::endl;
   _D.print("Underlying dictionary");
};

template <typename T>
void ConvolveDictionary<T>::multTrans(const Vector<T>& x, Vector<T>& b,
      const T alpha, const T beta) const {
   INTM size_feature_map=_sizeMapX*_sizeMapY;
   b.resize(_p);
   Image<T> I(x.rawX(),_w,_h,_nim);
   Matrix<T> patches;
   I.extractPatches(patches,_sizeEdges);
   Matrix<T> B(b.rawX(),_D.n(),size_feature_map);
   _D.mult(patches,B,true,false,alpha,beta);
}

template <typename T>
void ConvolveDictionary<T>::mult(const Vector<T>& x, Vector<T>& b,
      const T alpha, const T beta) const {
   INTM size_feature_map=_sizeMapX*_sizeMapY;
   b.resize(_w*_h*_nim);
   Matrix<T> tmp(x.rawX(),_D.n(),size_feature_map);
   Matrix<T> patches(_sizeEdges*_sizeEdges*_nim,size_feature_map);
   patches.setZeros();
   _D.mult(tmp,patches,false,false,alpha);
   Image<T> I(b.rawX(),_w,_h,_nim);
   I.combinePatches(patches,beta,1,false);
}

template <typename T>
void ConvolveDictionary<T>::mult(const SpVector<T>& x, Vector<T>& b, 
      const T alpha, const T beta) const {
   int size_feature_map=_sizeMapX*_sizeMapY;
   b.resize(_w*_h*_nim);
   SpMatrix<T> spalpha;
   x.toSpMatrix(spalpha,_D.n(),size_feature_map);
   Matrix<T> patches(_sizeEdges*_sizeEdges*_nim,size_feature_map);
   patches.setZeros();
   _D.mult(spalpha,patches,false,false,alpha);
   Image<T> I(b.rawX(),_w,_h,_nim);
   I.combinePatches(patches,beta,1,false);
}


#endif
