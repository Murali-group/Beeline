% 
% Usage:   alpha =mexConvFistaFlat(I,D,alpha0,param);
%
% Name: mexConvFistaFlat
%
% Description: performs convolutional sparse encoding of an
%              image I with a local dictionary D, using FISTA
%              and similar options as mexFistaFlat
%
% Inputs: I:  double nx x ny x nchannels   
%         D:  dictionary  
%         alpha0: initial weights  
%
% Output: alpha: output coefficients
%
% Author: Julien Mairal, 2014


