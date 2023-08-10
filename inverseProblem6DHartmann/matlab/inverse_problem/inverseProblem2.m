%%% Finding displaced log-normal pdf with maximum in 'Xmax' 
%%% and a tail such that P(x>'Xtail')='tailWeight'.
%%% The evaluation in values less than 'Xzero' is 0.
%%% Then the code generates 'nSample' samples of this distribution.
%%%

global Xmax;
global tailWeight;
global Xzero;
global Xtail;
nSample    = 1e3;
Xmax       = -3.322;
tailWeight = 0.1;
Xzero      = Xmax - 0.1;
Xtail      = -3;

x = [0; 0.1];
while norm(f(x)) > 1e-12
    x = x - Jacobian(x) \ (f(x).');
end

mu    = x(1);
sigma = x(2);

xPlot = Xzero+1e-14:1e-2:1;
%yPlot = arrayfun(@(x) logNormal(x,mu,sigma,Xzero), xPlot);
yPlot = lognpdf(xPlot-Xzero,mu,sigma);
%
sample = random('Lognormal',mu,sigma,[1 nSample]) + Xzero;
%edges  = [Xzero-0.2:1e-3:2]
%hist   = histogram(sample,edges);
%
%plot(xPlot,yPlot);
%hist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameters posterior sampling %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% metropolis-hastings sample %%
sigma_model = 1e-2;
covmatrixProposal = 1e-2*eye(6);
logpdf     = @(p) mahalanobis(h(p),Xmax-0.2,sigma_model);
proppdf    = @(p,q) mvnpdf(p,q,covmatrixProposal); 
proprnd    = @(p) rand([1 6]);
nsamplesMH = 1e5;
startPoint = 0.5*ones([1 6]);
postPredSample = mhsample(startPoint,nsamplesMH,'logpdf',logpdf,'proprnd',proprnd,'proppdf',proppdf);
sampleOverHartmann = zeros([1 nsamplesMH]);
for i = 1:nsamplesMH
  sampleOverHartmann(i) = h(postPredSample(i,:));
end

%%%%%%%%%%%%%%%
%%%% plots %%%%
%%%%%%%%%%%%%%%
edges  = [Xzero-0.2:0.05:0];
% posterior predictive
%display(histogram(sampleOverHartmann,edges,'Normalization','pdf'));
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
% prior
histogram(samplePriorOverHartmann,edges,'Normalization','pdf','FaceColor','blue'); hold on;
% target
histogram(sample,edges,'Normalization','pdf','FaceColor','green'); hold on;
% posterior predictive
histogram(sampleOverHartmann,edges,'Normalization','pdf','FaceColor','red');

