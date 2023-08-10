clear all; close all; clc;
addpath('../posteriorPredictiveApproachLeastSquare','-end')
addpath('../target_pdf','-end')

%nmin = 1;
%nmax = 30
%errors  = zeros(1,nmax-nmin+1);
%nparam  = zeros(1,nmax-nmin+1);
%for npar = nmin:nmax
%  nparam(npar-nmin+1) = npar;

% Number of legendre polynomials
npar = 17;
%% MODEL
% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
% Model approach functions
modelApproachFunctions = @(x) legendreApproachFunctionsRecursive(x,npar);
%% POSTERIOR PREDICTIVE APPROACH
[postpred,err] = posteriorPredictiveApproachLeastSquare(Xzero, 2, 1000, tpdf, modelApproachFunctions, true);
  
%  display(sprintf('%d -> %e',npar,err));
%  errors(npar-nmin+1) = err; 
%
%end

%plot(nparam, errors)
