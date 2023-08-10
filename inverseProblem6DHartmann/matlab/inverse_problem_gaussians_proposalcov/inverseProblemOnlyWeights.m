clear all; close all; clc;
addpath('posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
%% PARAMETERS
npar = 6;
for i=1:npar
  %            'name'                     initial      interval
  params{i} = {sprintf('\\theta_{%d}',i),   0.5    ,    0 , 10};
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
m1 = linspace(Xzero,-2,round(8*npar/10));
m2 = linspace(-2,2,npar-round(8*npar/10)+1);
m = [m1, m2(2:end)];
l = 1e-1*ones(1,npar);
modelApproachFunEval = gaussianFunctions(xres,m,l).';
modelerr = @(theta) sum((normalizedTrueDensity - theta*modelApproachFunEval).^2);
modelSigma2            = 0.01^2;
%% OPTIONS
options.method      = 'am';
options.nsimu       = 100000;
options.updatesigma = 0;
nord                = length(params);
% Proposal covariance
options.qcov        = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
burning             = 1000;
thin                = 10;
%% POSTERIOR PREDICTIVE APPROACH
[res0, chain, s2chain] = posteriorPredictiveApproach(xres, modelerr, modelSigma2, params, options, burning, thin, true);
for ii=1:size(chain,1)
  errchain(ii) = modelerr(chain(ii,:));
end
figure()
plot(errchain);
meanPredictedDensity = mean(chain,1)*modelApproachFunEval;
normalizedMeanPredictedDensity = meanPredictedDensity./sum(meanPredictedDensity);
lastPredictedDensity = chain(end,:)*modelApproachFunEval;
normalizedLastPredictedDensity = lastPredictedDensity./sum(lastPredictedDensity);
figure()
plot(xres,normalizedMeanPredictedDensity, '-r'); hold on;
plot(xres,normalizedLastPredictedDensity, '-g');
plot(xres,normalizedTrueDensity, '-b')
hold off

