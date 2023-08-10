clear all; close all; clc;
addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
%% PARAMETERS
npar = 17;
for i=1:npar
  %            'name'                     initial      interval
  params{i} = {sprintf('\\theta_{%d}',i),   0.5    ,  -10 , 10};
end
%% MODEL
% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
% Points
nres = 1000;
xres = linspace(Xzero,2,nres);
% Normalized target density evaluations
tpdfOVERxres          = tpdf(xres);
normalizedTrueDensity = tpdfOVERxres./sum(tpdfOVERxres);
% Model approach functions evaluations
modelApproachFunctions = @(x) legendreApproachFunctionsRecursive(x,npar);
modelSigma2            = 0.01^2;
%% OPTIONS
options.method      = 'am';
options.nsimu       = 10000;
options.updatesigma = 0;
nord                = length(params);
% Proposal covariance
options.qcov        = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
burning             = 1000;
thin                = 10;
%% POSTERIOR PREDICTIVE APPROACH
predictedTheta = posteriorPredictiveApproach(xres, normalizedTrueDensity, modelApproachFunctions, modelSigma2, params, options, burning, thin, true, false);
