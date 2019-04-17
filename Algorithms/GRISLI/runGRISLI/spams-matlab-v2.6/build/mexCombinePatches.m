% 
% Usage:   I =mexCombinePatches(X,I0,n,lambda,step,normalize);
%
% Name: mexCombinePatches
%
% Description: Combine patches extracted with mexExtractPatches into a new image 
%              I = lambda I0 + combined(X) if normalize, averaging is
%              performed; otherwise the patches are just summed.
%
% Inputs: I0:  double nx x ny x nchannels   
%
% Output: I: double nx x ny x nchannels 
%
% Author: Julien Mairal, 2014


