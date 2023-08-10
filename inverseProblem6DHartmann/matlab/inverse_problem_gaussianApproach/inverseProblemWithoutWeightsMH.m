clear all; close all; clc;
addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
%% TARGET DENSITY
% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
%% MODEL
% Points
nres = 1000;
xres = linspace(Xzero,2,nres);
% Normalized target density evaluations
tpdfOVERxres          = tpdf(xres);
normalizedTrueDensity = tpdfOVERxres./sum(tpdfOVERxres);
% Model approach functions evaluations
ngaussians = 20;
modelApproachFunctions = @(theta,x) gaussianApproachFunctionsWithoutWeightsMH(theta,2*ngaussians,x);
modelSigma2            = 0.01^2;
%% PARAMETERS
%f = @(theta) sum((normalizedTrueDensity - gaussianApproachFunctionsMH(theta,3*ngaussians,xres)).^2);
%for i=1:ngaussians
%  theta0(i) = 1/ngaussians;
%end
%for i=1:ngaussians
%  theta0(ngaussians+i) = Xzero + i*(10/ngaussians);
%end
%for i=1:ngaussians
%  theta0(2*ngaussians+i) = 1e-2;
%end
%theta  = fminsearch(f,theta0);
%theta = [ -3, -3, -3, -3, 1e-2, 1e-2, 1e-2, 1e-2 ];
for i=1:ngaussians
  %                         'name'            initial                    interval
  params{i} = {sprintf('\\thetaMean_{%d}',i), -4+(i/ngaussians)*6 ,  -5 , 3};
end
for i=1:ngaussians
  %                                    'name'               initial    interval
  params{ngaussians+i} = {sprintf('\\thetaVariance_{%d}',i), 1e-2 ,  1e-10 , 1};
end
%% OPTIONS
options.method      = 'am';
options.nsimu       = 100000;
options.updatesigma = 0;
nord                = length(params);
% Proposal covariance
%options.qcov        = eye(nord).*0.1^2;
options.qcov = zeros(nord);
%options.qcov        = (ones(nord)-eye(nord)).*0.01.^2;
for i=1:ngaussians
  % means variances
  options.qcov(i,i) = 0.01;
end
for i=1:ngaussians
  % variances variances
  options.qcov(ngaussians+i,ngaussians+i) = 1e-3;
end
burning             = 10000;
thin                = 10;
%% POSTERIOR PREDICTIVE APPROACH
[chain, predictedTheta] = posteriorPredictiveApproachMH(xres, normalizedTrueDensity, modelApproachFunctions, modelSigma2, params, options, burning, thin, true, false);

%function a = f2(x,m,l,tdensity)  
%aux = x.'-m;
%H   = (aux.*aux)./(2*l);
%a   = H \ tdensity;
%end
