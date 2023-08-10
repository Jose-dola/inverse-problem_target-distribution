clear all; close all; clc;
addpath('posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
% Number of gaussians
ngaussians = 4;
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
modelerr = @(theta) errGaussians(xres, ...
                       normalizedTrueDensity, ... 
                       [1 theta(1:ngaussians-1)], ...
                       theta(ngaussians:2*ngaussians-1), ...
                       theta(2*ngaussians:end)); 
modelSigma2 = 0.01^2;
%% PARAMETERS
m = linspace(Xzero,2,ngaussians);
l = 1e-1*ones(1,ngaussians);
bestTheta = [1 fminsearch(modelerr,[ ones(1,ngaussians-1) m l ] )];
theta0 = [ ones(1,ngaussians) m l ];
for i=1:ngaussians-1
  %            'name'                              initial            interval
  params{i} = {sprintf('\\thetaWeight_{%d}',i),   theta0(i+1)    ,   -10 , 10};
end
for i=1:ngaussians
  %                   'name'                                     initial                interval
  params{ngaussians-1+i} = {sprintf('\\thetaMean_{%d}',i), theta0(ngaussians+i)    ,    Xzero, 2 };
end
for i=1:ngaussians
  %                   'name'                                     initial                  interval
  params{2*ngaussians-1+i} = {sprintf('\\thetaVar_{%d}',i), theta0(2*ngaussians+i)    ,  1e-10, 1 };
end
%% OPTIONS
options.method      = 'am';
options.nsimu       = 3*1e4;
options.updatesigma = 0;
nord                = length(params);
% Proposal covariance
options.qcov        = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
burning             = 1e3;
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
bestChain = [1 chain(j,:)];

% plot chain error
figure()
plot(errchain);

% plot true density and predicted densities for different theta values
figure()

meanChain = [1 mean(chain,1)];
normalizedMeanPredictedDensity = predictedDensity(xres,meanChain,ngaussians);
plot(xres,normalizedMeanPredictedDensity, '-r'); 

hold on;

lastChain = [1 chain(end,:)];
normalizedLastPredictedDensity = predictedDensity(xres,lastChain,ngaussians);
plot(xres,normalizedLastPredictedDensity, '-g');

normalizedTheta0PredictedDensity = predictedDensity(xres,theta0,ngaussians);
plot(xres,normalizedTheta0PredictedDensity, '-c');

normalizedBestThetaPredictedDensity = predictedDensity(xres,bestTheta,ngaussians);
plot(xres,normalizedBestThetaPredictedDensity, '-m');

normalizedBestChainPredictedDensity = predictedDensity(xres,bestChain,ngaussians);
plot(xres,normalizedBestChainPredictedDensity, '-k');

plot(xres,normalizedTrueDensity, '-b');

hold off;

function e = errGaussians(x,normalizedTrueDensity,alpha,m,l)
  predDensity = alpha*gaussianFunctions(x,m,l).';
  e = sum((normalizedTrueDensity - predDensity./sum(predDensity)).^2);
  %e = sum(abs(normalizedTrueDensity - predDensity./sum(predDensity)));
end

function pd = predictedDensity(x,theta,ngaussians)
  alpha = theta( 1             : ngaussians   );
  m     = theta( ngaussians+1  : 2*ngaussians );
  l     = theta( 2*ngaussians+1: end          );
  pd = alpha*gaussianFunctions(x,m,l).';
  pd = pd./sum(pd); 
end
