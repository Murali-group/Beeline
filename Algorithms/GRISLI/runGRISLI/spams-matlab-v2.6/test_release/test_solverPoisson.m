rand('seed',0);
randn('seed',0);
format compact;
randn('seed',0);
param.numThreads=-1; % all cores (-1 by default)
param.verbose=true;   % verbosity, false by default
param.lambda=0.05; % regularization parameter
param.it0=10;      % frequency for duality gap computations
param.max_it=200; % maximum number of iterations
param.L0=0.01;
param.tol=1e-6;
param.delta=1e-5;

X=rand(1000,2000);
X=mexNormalize(X);
Y=rand(1000,1);
Y=mexNormalize(Y);
W0=zeros(size(X,2),size(Y,2));
% Regression experiments 
% 100 regression problems with the same design matrix X.
fprintf('\nVarious regression experiments\n');
fprintf('\nFISTA + Regression l1\n');
param.regul='l1';
tic
[W optim_info]=solverPoisson(Y,X,W0,param);
t=toc;
fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));


