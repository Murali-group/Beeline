clear all;
mex -v -DNDEBUG -I./linalg/ -I./decomp/ -lgomp 'decomp/mex/mexLasso.cpp' -o 'build/mexLasso.mex';
