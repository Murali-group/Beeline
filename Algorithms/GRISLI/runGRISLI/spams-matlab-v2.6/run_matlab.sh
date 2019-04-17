#!/bin/sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/
matlab $* -singleCompThread -r "addpath('./build/'); addpath('./test_release'); setenv('MKL_NUM_THREADS','1'); setenv('MKL_SERIAL','YES');"
