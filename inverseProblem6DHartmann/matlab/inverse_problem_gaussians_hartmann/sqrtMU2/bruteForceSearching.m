clear all; close all; clc;

addpath('../posteriorPredictiveApproach','-end')
addpath('../../target_pdf','-end')
% Number of gaussians
ngaussians = 3;
% Parameters Dimension
parametersDimension = 6;
%% Target Density
Xmax       = -3.322;
Xzero      = Xmax - 0.1;
Xtail      = -2;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
%% Histo Structure
h = histo(Xzero,2,100);
%% true density
trueDensity = tpdf( h.intervalrepresentativePoints() );
%% mcmc model error
modelerr = @(theta) E(theta,ngaussians,parametersDimension,h,1e3,trueDensity);
modelSigma2 = 0.01^2;


minerr = 1e10;
iter   = 1;
while(true)

display(sprintf('\n\n============')); 
display(sprintf('%d ITERATION',iter)); 
display(sprintf('============\n'));

clear options;
clear params; 
clear burning; 
clear thin;
clear nord;
clear proposalCov;
clear errc;
clear c;
clear s2c;
clear r0;
clear t0;
clear theta0;
clear res0;
clear chain; 
clear s2chain;
clear bestchainerr;
clear errchain;
clear fi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finding the parameters that result in Xzero when we apply the hartmann function
modelerr2    = @(x) (hart6(x)-Xmax)^2;
modelSigma2 = 0.01^2;
for i =1:6
  par{i} = {sprintf('\\theta_{%d}',i),   rand()    ,   0 , 1};
end
options.method      = 'am';
options.nsimu       = 2*1e4;
options.updatesigma = 0;
nord                = length(par);
% Proposal covariance
options.qcov        = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
burning             = options.nsimu-round(options.nsimu/5);
thin                = 1;
[r0, c, s2c] = posteriorPredictiveApproach(modelerr2, modelSigma2, par, options, burning, thin, false);
j = 1;
errc(1) = modelerr2(c(1,:));
for i=2:size(c,1)
  errc(i) = modelerr2(c(i,:));
  if (errc(i) < errc(j)) j = i; end
end
t0 = c(j,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
nord = 2*ngaussians+parametersDimension*ngaussians-1;
proposalCov = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
theta0(1)=1/2;
theta0(2)=1/2;
for i=1:ngaussians-1
  %                           'name'               initial        interval
  params{i} = {sprintf('\\thetaWeight_{%d}',i),     theta0(i)  ,   0 , 1};
  proposalCov(i,i) = 100;
end
for i = ngaussians : ngaussians+parametersDimension*ngaussians-1
  n = floor((i-ngaussians)/parametersDimension)+1;
  index = mod((i-ngaussians),parametersDimension)+1;
  if(index  == 1) j=1; end
  theta0(i) = t0(j);
  %                             'name'                 initial        interval
  params{i} = {sprintf('\\thetaMean_{%d,%d}',n,index), theta0(i) ,    Xzero, 2 };
  proposalCov(i,i) = 100;
  j = j+1;
end
for i = ngaussians+parametersDimension*ngaussians : 2*ngaussians+parametersDimension*ngaussians-1
  j = i-(ngaussians+parametersDimension*ngaussians)+1;
  theta0(i) = rand()*1e-3;
  %                   'name'                  initial         interval
  params{i} = {sprintf('\\thetaVar_{%d}',j),  theta0(i) ,  1e-10, 1 };
  proposalCov(i,i) = 0.05^2;
end
%% OPTIONS
options.method      = 'am';
options.nsimu       = 800;
options.updatesigma = 0;
% Proposal covariance
options.qcov        = proposalCov;
burning             = 1;
thin                = 1;
%% POSTERIOR PREDICTIVE APPROACH
[res0, chain, s2chain] = posteriorPredictiveApproach(modelerr, modelSigma2, params, options, burning, thin, false);

% errors and best approximation in the chain
k = options.nsimu-round(0.9*options.nsimu);
bestChainIndex = 1;
bestchainerr = 1e10;
for i=k:size(chain,1)
  errchain = modelerr(chain(i,:));
  if (errchain < bestchainerr) 
    bestchainerr = errchain; 
    bestChainIndex = i;
  end
end

if (bestchainerr < minerr)
  minerr = bestchainerr;
  % saved the chain
  fi = fopen( sprintf('/home/josedorado/data/chain.doubleMatrix800x23_%d.bin',iter) , 'w' );
  fwrite(fi, chain, 'double');   
  fclose(fi); 
  fi = fopen( sprintf('/home/josedorado/data/theta0.doubleVector1x23_%d.bin',iter) , 'w' );
  fwrite(fi, theta0, 'double');   
  fclose(fi); 
  fi = fopen( sprintf('/home/josedorado/data/bestChain.doubleVector1x23_%d.bin',iter) , 'w' );
  fwrite(fi, chain(bestChainIndex,:), 'double');   
  fclose(fi); 
end

iter = iter+1;
end


