%%% Finding displaced log-normal pdf with maximum in 'Xmax' 
%%% and a tail such that P(x>'Xtail')='tailWeight'.
%%% The evaluation in values less than 'Xzero' is 0.
%%% Then the code generates 'nSample' samples of this distribution.
%%%

global Xmax;
global tailWeight;
global Xzero;
global Xtail;
nSample    = 1e6;
Xmax       = -3.3592;
tailWeight = 0.1;
Xzero      = Xmax - 0.1;
Xtail      = -2;

x = [0; 0.1];
while norm(f(x)) > 1e-12
    x = x - Jacobian(x) \ (f(x).');
end

mu    = x(1);
sigma = x(2);

%xPlot = Xzero+1e-14:1e-2:1;
%yPlot = arrayfun(@(x) logNormal(x,mu,sigma,Xzero), xPlot);
%yPlot = lognpdf(xPlot-Xzero,mu,sigma);
%
sample              = random('Lognormal',mu,sigma,[1 nSample]) + Xzero;

%edges  = [Xzero-0.2:1e-3:2]
%hist   = histogram(sample,edges);
%
%plot(xPlot,yPlot);
%hist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameters' posterior sampling %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% metropolis-hastings sample (MCMC toolbox) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../mcmcstat');
clear data model options params

%% data structure %%
data.sumOfTerms          = sum(sample);
data.sumOfQuadraticTerms = sum(sample.*sample);
data.nTerms              = nSample;
%% parameters %%
startPoint = [0.202 0.150 0.477 0.275 0.312 0.657];
%startPoint = 0.5*ones([1 6]);
for i=1:6
  params{i} = {sprintf('theta_%d',i),startPoint(i),0,1};
end
%% model %%
%% -2*log(likelihood) -> sum of squares
data.sigma_model = 1e-1;
model.ssfun=@(x,data) sumOfSquares(h(x),data);
%% number of observations
model.N = nSample;
%% initial error variance 
%model.sigma2 = 1e-3;
%% options %%
%% number of simulations
options.nsimu = 5*1e3;
options.method = 'am';
%% update sigma2
options.updatesigma = 0;
%% proposal covariance
options.qcov = 1e-2.*eye(6);
%% get samples from mcmc toolbox
[res,chain] = mcmcrun(model,data,params,options);

sampleMHOverHartmann = zeros([1 length(chain)]);
for i = 1:length(chain)
  sampleMHOverHartmann(i) = h(chain(i,:));
end

%%% metropolis-hastings sample %%
%sigma_model = 1e-3;
%%covmatrixProposal = 1e-3*eye(6);
%sigma_proposal = 1e-3;
%%logpdf     = @(p)   mahalanobisDistance(h(p),sample,sigma_model);
%logpdf     = @(p)   mahalanobisDistanceToSample(h(p),nSample,sumOfTerms, sumOfQuadraticTerms, sigma_model);
%%proppdf    = @(p,q) mvnpdf(p,q,covmatrixProposal); 
%proppdf    = @(p,q) 1; 
%proprnd    = @(p)   proposal(p,sigma_proposal);
%nsamplesMH = 1e3;
%%startPoint = [0.202 0.150 0.477 0.275 0.312 0.657];
%startPoint = 0.5*ones([1 6]);
%%h(startPoint)
%postPredSample = mhsample(startPoint,nsamplesMH,'logpdf',logpdf,'proprnd',proprnd,'proppdf',proppdf);
%sampleMHOverHartmann = zeros([1 nsamplesMH]);
%for i = 1:nsamplesMH
%  sampleMHOverHartmann(i) = h(postPredSample(i,:));
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plots and save the results %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges  = [Xzero-0.2:0.05:0.5];
% posterior predictive
%display(histogram(sampleMHOverHartmann,edges,'Normalization','pdf'));
% target
%display(histogram(sample,edges));
% prior
nsampleprior = 1e5;
samplePrior = rand([nsampleprior 6]);
samplePriorOverHartmann = zeros([1 nsampleprior]);
for i = 1:nsampleprior
  samplePriorOverHartmann(i) = h(samplePrior(i,:));
end
figure()
%% prior
fi = fopen('samplePriorOverHartmann.bin','w');
fwrite(fi, samplePriorOverHartmann);   
fclose(fi); 
histogram(samplePriorOverHartmann,edges,'Normalization','pdf','FaceColor','blue'); hold on;
%% target
fi = fopen('sample.bin','w');
fwrite(fi, sample);   
fclose(fi); 
histogram(sample,edges,'Normalization','pdf','FaceColor','green'); hold on;
%% posterior predictive
fi = fopen('sampleMHOverHartmann.bin','w');
fwrite(fi, sampleMHOverHartmann);   
fclose(fi); 
histogram(sampleMHOverHartmann,edges,'Normalization','pdf','FaceColor','red');

figure()
plot(chain)

