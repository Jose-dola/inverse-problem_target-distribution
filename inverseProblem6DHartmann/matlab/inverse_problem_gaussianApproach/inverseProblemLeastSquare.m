clear all; close all; clc;
addpath('../posteriorPredictiveApproachLeastSquare','-end')
addpath('../target_pdf','-end')

%% MODEL
% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
% Model approach functions
% means' vector
nmeans = 100;
m = linspace(Xzero,2,nmeans);
% variances' vector
nvariances = 100;
l = linspace(1e-3,1,nvariances);
%nvariances = 1;
%l = linspace(1-0.1,1,nvariances);
l = l.*l;
modelApproachFunctions = @(x) gaussianApproachFunctions(x,m,l);
%% POSTERIOR PREDICTIVE APPROACH
[postpred,err] = posteriorPredictiveApproachLeastSquare(Xzero, 2, 1e4, tpdf, modelApproachFunctions, true);
err 
