clear all; close all; clc;
addpath('posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
% Number of gaussians
ngaussians = 5;
% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
% Points
nres = 1000;
xres = linspace(Xzero,2,nres);
%% MODEL
% Normalized target density evaluations
tpdfOVERxres          = tpdf(xres);
normalizedTrueDensity = tpdfOVERxres./sum(tpdfOVERxres);
% Model approach functions evaluations
l = 1e-3*ones(1,ngaussians);
modelerr = @(theta) errGaussians(xres,normalizedTrueDensity,[1 theta(1:ngaussians-1)],theta(ngaussians:end),l); 
modelSigma2 = 0.01^2;
%% PARAMETERS
%m1 = linspace(Xzero,-2,round(8*ngaussians/10));
%m2 = linspace(-2,2,ngaussians-round(8*ngaussians/10)+1);
%m = [m1, m2(2:end)];
m = linspace(Xzero,2,ngaussians);
bestTheta = [1 fminsearch(modelerr,[ ones(1,ngaussians-1) m ] )];
theta0 = [ ones(1,ngaussians) m ];
%gradient = @(theta) gradientModelerr(theta, l, gaussianFunctions(xres,m,l), ngaussians, xres, normalizedTrueDensity);
%theta0 = [ (1/ngaussians)*ones(1,ngaussians) m ];
%theta0 = gradientDescent(theta0, modelerr, gradient, 1e-3, 1e100, 1e-5);
for i=1:ngaussians-1
  %            'name'                              initial            interval
  params{i} = {sprintf('\\thetaWeight_{%d}',i),   theta0(i+1)    ,           0 , 10};
end
for i=1:ngaussians
  %                   'name'                                     initial                  interval
  params{ngaussians-1+i} = {sprintf('\\thetaMean_{%d}',i), theta0(ngaussians+i)    ,    Xzero, 2 };
end
%% OPTIONS
options.method      = 'am';
options.nsimu       = 2*1e4;
options.updatesigma = 0;
nord                = length(params);
% Proposal covariance

burning             = options.nsimu-round(options.nsimu/5);
thin                = 1;
%% POSTERIOR PREDICTIVE APPROACH
[res0, chain, s2chain] = posteriorPredictiveApproach(xres, modelerr, modelSigma2, params, options, burning, thin, true);
%% RESULTS %%
% errors and best approximation in the chain
j = 1;
errchain(1) = modelerr(chain(1,:));
for i=2:size(chain,1)
  errchain(i) = modelerr(chain(i,:));
  if (errchain(i) < errchain(j)) j = i; end
end

% plot chain error
figure()
plot(errchain);

% plot true density and predicted densities for different theta values
figure()

meanChain = [1 mean(chain,1)];
meanPredictedDensity = meanChain(1:ngaussians)*gaussianFunctions(xres,meanChain(ngaussians+1:end),l).';
normalizedMeanPredictedDensity = meanPredictedDensity./sum(meanPredictedDensity);
plot(xres,normalizedMeanPredictedDensity, '-r'); 

hold on;

lastChain = [1 chain(end,:)];
lastPredictedDensity = lastChain(1:ngaussians)*gaussianFunctions(xres,lastChain(ngaussians+1:end),l).';
normalizedLastPredictedDensity = lastPredictedDensity./sum(lastPredictedDensity);
plot(xres,normalizedLastPredictedDensity, '-g');

theta0PredictedDensity = theta0(1:ngaussians)*gaussianFunctions(xres,theta0(ngaussians+1:end),l).';
normalizedTheta0PredictedDensity = theta0PredictedDensity./sum(theta0PredictedDensity);
plot(xres,normalizedTheta0PredictedDensity, '-c');

bestThetaPredictedDensity = bestTheta(1:ngaussians)*gaussianFunctions(xres,bestTheta(ngaussians+1:end),l).';
normalizedBestThetaPredictedDensity = bestThetaPredictedDensity./sum(bestThetaPredictedDensity);
plot(xres,normalizedBestThetaPredictedDensity, '-m');

bestChain = [1 chain(j,:)];
bestChainPredictedDensity = bestChain(1:ngaussians)*gaussianFunctions(xres,bestChain(ngaussians+1:end),l).';
normalizedBestChainPredictedDensity = bestChainPredictedDensity./sum(bestChainPredictedDensity);
plot(xres,normalizedBestChainPredictedDensity, '-k');

plot(xres,normalizedTrueDensity, '-b');

hold off;

function e = errGaussians(x,normalizedTrueDensity,alpha,m,l)
  predDensity = alpha*gaussianFunctions(x,m,l).';
  e = sum((normalizedTrueDensity - predDensity./sum(predDensity)).^2);
  %e = sum(abs(normalizedTrueDensity - predDensity./sum(predDensity)));
end

function g = gradientModelerr(theta, l, exps, ngaussians, points, td)
  alpha = theta(1:ngaussians);
  m     = theta(ngaussians+1 : end);
%  exps  = gaussianFunctions(points,m,l);
  aux1  = td - alpha*exps.';
  aux2  = (points.'-m)./l;
  g = [ 2*aux1*exps   2*points*(exps.*aux2.*alpha) ];
end
