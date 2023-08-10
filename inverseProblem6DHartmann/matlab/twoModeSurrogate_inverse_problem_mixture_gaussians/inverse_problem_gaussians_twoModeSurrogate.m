clear all; close all; clc;

%numberOfGaussians = 3;
%parametersDimension = 6;
%theta = rand(1,numberOfGaussians-1+numberOfGaussians*parametersDimension+numberOfGaussians);
%for i = 0:numberOfGaussians-1
%  theta(end-i) = theta(end-i)*1e-2;
%end
%[weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension)
%h = histo(-5,5,20);
%hUniform = histo(-5,5,20);
%for i = 1:1e3
%  h.pointToHisto( getRandomFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension) );
%  hUniform.pointToHisto( rand(1,6) );
%end
%h.normalize();
%hUniform.normalize();
%figure()
%h.plot(); hold on;
%hUniform.plot();
%hold off;

%td = normpdf(h.intervalrepresentativePoints(),0,1);
%L2error(h,td)

addpath('posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
% Number of gaussians
ngaussians = 3;
% Parameters Dimension
parametersDimension = 2;
%% Target Density
Xmax       = -1.9;
Xzero      = Xmax - 0.2;
Xtail      = -1;
tailWeight = 0.1;
[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histo Structure
h = histo(Xzero,1,100);
%% true density
trueDensity = tpdf( h.intervalrepresentativePoints() );
%% mcmc model error
modelerr = @(theta) E(theta,ngaussians,parametersDimension,h,1e3,trueDensity);
modelSigma2 = 0.01^2;
%% Parameters
%theta0 = rand(1,ngaussians-1+ngaussians*parametersDimension+ngaussians);
%for i = 0:ngaussians-1
%  theta0(end-i) = theta0(end-i)*1e-2;
%end
nord = 2*ngaussians+parametersDimension*ngaussians-1;
proposalCov = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
theta0(1)=1/2;
theta0(2)=1/2;
for i=1:ngaussians-1
  %theta0(i) = 0.5;
  %                           'name'               initial        interval
  params{i} = {sprintf('\\thetaWeight_{%d}',i),     theta0(i)  ,   0 , 1};
  proposalCov(i,i) = 0.05^2;
end
for i = ngaussians : ngaussians+parametersDimension*ngaussians-1
  n = floor((i-ngaussians)/parametersDimension)+1;
  index = mod((i-ngaussians),parametersDimension)+1;
  if(index  == 1) j=1; end
  theta0(i) = rand;
  %                             'name'                 initial        interval
  params{i} = {sprintf('\\thetaMean_{%d,%d}',n,index), theta0(i) ,    Xzero, 2 };
  proposalCov(i,i) = 0.005;
  j = j+1;
end
for i = ngaussians+parametersDimension*ngaussians : 2*ngaussians+parametersDimension*ngaussians-1
  j = i-(ngaussians+parametersDimension*ngaussians)+1;
  theta0(i) = rand()*5*1e-2;
  %                   'name'                  initial         interval
  params{i} = {sprintf('\\thetaVar_{%d}',j),  theta0(i) ,  0, 1 };
  proposalCov(i,i) = 0.05^2;
end
%% OPTIONS
options.method      = 'am';
options.nsimu       = 1000;
options.updatesigma = 0;
%nord                = length(params);
% Proposal covariance
% options.qcov        = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
options.qcov        = proposalCov;
%burning             = options.nsimu-round(0.7*options.nsimu);
burning             = 1;
thin                = 1;
%% POSTERIOR PREDICTIVE APPROACH
[res0, chain, s2chain] = posteriorPredictiveApproach(modelerr, modelSigma2, params, options, burning, thin, false);

%% RESULTS %%
% saved the chain
fi = fopen('chain.doubleMatrix.bin','w');
fwrite(fi, chain, 'double');   
fclose(fi); 
fi = fopen('chainSize.intVector.bin','w');
fwrite(fi, size(chain), 'double');   
fclose(fi); 

% errors and best approximation in the chain
j = 1;
errchain(1) = modelerr(chain(1,:));
for i=2:size(chain,1)
  errchain(i) = modelerr(chain(i,:));
  if (errchain(i) < errchain(j)) j = i; end
end
bestChain = chain(j,:);
fi = fopen('bestChainElemntError.double.bin','w');
fwrite(fi, modelerr(bestChain), 'double'); 
fclose(fi);  

% plot chain error
figure()
plot(errchain);

% plot true density and predicted densities for different theta values
figure()

meanChain = mean(chain,1);
fillHisto(meanChain,ngaussians,parametersDimension,h,1e3);
h.plot('-r');

hold on;

lastChain = chain(end,:);
fillHisto(lastChain,ngaussians,parametersDimension,h,1e3);
h.plot('-g');

fillHisto(theta0,ngaussians,parametersDimension,h,1e3);
h.plot('-c');

fillHisto(bestChain,ngaussians,parametersDimension,h,1e3);
h.plot('-k');

plot(h.intervalrepresentativePoints,trueDensity, '-b');

hold off;

% plot how the means move
figure()
p = 0:0.005:1;
[X Y Z] = twoModeSurrogateGrid(p,p);
pcolor(X,Y,Z);             
shading flat;              
hold on;                   
plot(chain(:,3),chain(:,4));
plot(chain(:,5),chain(:,6));
plot(chain(:,7),chain(:,8));
hold off;


