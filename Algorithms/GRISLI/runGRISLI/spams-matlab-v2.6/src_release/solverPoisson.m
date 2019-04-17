% 
% Usage:   [W, [optim]]=solverPoisson(Y,X,W0,param);
%
% Name: solverPoisson
%
% Description: solverPoisson solves the regularized Poisson regression
% problem for every column of Y:
%
%   min_w  \sum_i (x_i' w) + delta  - y_i log( x_i' w + delta) + lambda psi(w) (1)
%
% delta should be positive. The method solves (1) for decreasing values of delta,
% until it reaches the desired target value. The algorithm ISTA is used
% with Barzilai-Borwein steps.
%
% Inputs: Y:  double m x n matrix (non-negative values)
%         X:  double m x p matrix (non-negative values)
%         W0:  double p x n matrix (non-negative values)
%         param: struct
%            param.tol, param.max_it, param.verbose, param.L0,
%            param.it0, param.lambda as for mexFistaFlat
%            param.regul,  as for mexProximalFlat
%
% Output:
%         W: double p x n matrix (non-negative values)

function [W optim] = solverPoisson(Y,X,W0,param)
param.delta=param.delta;
param.lambda=param.lambda;
param.ista=true;
param.intercept=false;
param.pos=true;
param.L0=1e-5;
param.linesearch_mode=2;
param.loss='poisson';
param.verbose=false;
tabdelta=logspace(0,log10(param.delta),1-log10(param.delta));
for delta = tabdelta
   param2=param;
   param2.delta=delta;
   param2.verbose=abs(delta-tabdelta(end)) < 1e-10 && param.verbose;
   [W optim]=mexFistaFlat(Y,X,W0,param2);
   W0=W;
end

