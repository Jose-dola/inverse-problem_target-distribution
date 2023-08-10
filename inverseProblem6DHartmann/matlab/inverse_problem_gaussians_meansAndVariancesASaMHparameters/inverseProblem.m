clear all; close all; clc;
addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
%% PARAMETERS
ngaussians = 10;
m1 = linspace(Xzero,-2,round(8*ngaussians/10));
m2 = linspace(-2,2,ngaussians-round(8*ngaussians/10)+1);
m = [m1, m2(2:end)];
l = 1e-1*ones(1,ngaussians);
for i=1:ngaussians
  %            'name'                         initial interval     
  params{i} = {sprintf('\\thetaWheigh_{%d}',i),   1,  -10, 10};
end
for i=1:ngaussians
  %            'name'                                   initial  interval
  params{ngaussians+i} = {sprintf('\\thetaMean_{%d}',i), m(i) ,  -5 , 3};
end
%for i=1:ngaussians
  %                          'name'                           initial   interval
  %params{2*ngaussians+i} = {sprintf('\\thetaVariance_{%d}',i), 1e-2 ,  1e-10 , 1};
%end
%% MODEL
% Points
nres = 1000;
xres = linspace(Xzero,2,nres);
% Normalized target density evaluations
tpdfOVERxres          = tpdf(xres);
normalizedTrueDensity = tpdfOVERxres./sum(tpdfOVERxres);
%modelApproachFunctions = @(x) (2*(2-Xzero)/nres)*gaussianFunctions(x,m,l);
modelApproachFunctions = @(x) gaussianFunctions(x,m,l);
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
