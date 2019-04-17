clear all; 
rand('seed',0);
I=double(imread('data/lena.png'))/255;
% extract 8 x 8 patches
X=im2col(I,[8 8],'sliding');
per=randperm(size(X,2));
X=X(:,per(1:10000));
X=X-repmat(mean(X),[size(X,1) 1]);
nrm=sqrt(sum(X.^2));
me=median(nrm);
X=X ./ repmat(max(nrm,me),[size(X,1) 1]);

param.p=64;  % learns a dictionary with 64 elements
param.robust=false;
param.epsilon=1e-3;  % width for Huber loss
param.computeXtX=true;
param.stepsFISTA=0; 
param.stepsAS=10;
param.numThreads=-1;

%%%%%%%%%% FIRST EXPERIMENT %%%%%%%%%%%
tic
[Z A B] = mexArchetypalAnalysis(X,param);
t=toc;
fprintf('time of computation for Archetypal Analysis: %f\n',t);
R=norm(X-Z*A,'fro')^2;
fprintf('objective function: %f\n',R);

fprintf('Re-Evaluating cost function after re-updating A...\n');
paramdecomp.computeXtX=true;
paramdecomp.numThreads=-1;
A2=mexDecompSimplex(X,Z,param);
R=norm(X-Z*A2,'fro')^2;
fprintf('objective function: %f\n',R);

%%%%%%%%%% Second EXPERIMENT %%%%%%%%%%%
tic
[Z2 A2] = mexArchetypalAnalysis(X,param,Z);
t=toc;
fprintf('time of computation for Archetypal Analysis: %f\n',t);
fprintf('Evaluating cost function...\n');
R=norm(X-Z2*A2,'fro')^2;
fprintf('objective function: %f\n',R);

%%%%%%%%%% THIRD EXPERIMENT %%%%%%%%%%%
tic
param.robust=true;
[Z3] = mexArchetypalAnalysis(X,param);
t=toc;
fprintf('time of computation for Robust Archetypal Analysis: %f\n',t);



