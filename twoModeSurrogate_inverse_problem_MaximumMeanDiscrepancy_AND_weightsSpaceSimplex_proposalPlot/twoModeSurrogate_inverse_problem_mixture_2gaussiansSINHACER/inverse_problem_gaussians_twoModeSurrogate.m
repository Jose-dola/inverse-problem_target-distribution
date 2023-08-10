clear all; close all; clc;
fsize=28; fname='times';

addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
addpath('../histoStructure','-end')
% Number of gaussians
ngaussians = 2;
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
modelSigma2 = 0.05;
modelPriorSigma2S20 = modelSigma2;
modelPriorSigma2N0  = 10000;
%% Parameters
nord = 2*ngaussians+parametersDimension*ngaussians-1;
proposalCov = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
for i=1:ngaussians-1
  theta0(i) = 1/2;
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
  params{i} = {sprintf('\\thetaMean_{%d,%d}',n,index), theta0(i) ,    0, 1 };
  proposalCov(i,i) = 0.0005;
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
options.nsimu       = 50000;
options.updatesigma = 1;
options.qcov        = proposalCov;
burning             = 20000;
thin                = 1;
%% POSTERIOR PREDICTIVE APPROACH -- RUNING MCMC
[res0, chain, s2chain] = posteriorPredictiveApproach(modelerr, modelSigma2, modelPriorSigma2S20, modelPriorSigma2N0, params, options, burning, thin, false);

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
lineWidth = 2;

meanChain = mean(chain,1);
fillHisto(meanChain,ngaussians,parametersDimension,h,1e3);
[f,xi] = h.ksdens;
plot(xi,f,'-r','LineWidth',lineWidth);

hold on;

lastChain = chain(end,:);
fillHisto(lastChain,ngaussians,parametersDimension,h,1e3);
[f,xi] = h.ksdens;
plot(xi,f,'-g','LineWidth',lineWidth);

fillHisto(theta0,ngaussians,parametersDimension,h,1e3);
[f,xi] = h.ksdens;
plot(xi,f,'-c','LineWidth',lineWidth);

fillHisto(bestChain,ngaussians,parametersDimension,h,1e3);
[f,xi] = h.ksdens;
plot(xi,f,'-k','LineWidth',lineWidth);

plot(h.intervalrepresentativePoints,trueDensity, '-b','LineWidth',lineWidth);

xlabel('y','FontName',fname,'fontsize',fsize);
ylabel('P(y)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');

hold off;

% plot how the means move

figure()
for i = 1:length(chain)
  weights(i,:)   = mcmcParamToRealParam(chain(i,:),ngaussians,parametersDimension);
  evalChain(i,1) = twoModeSurrogate([chain(i,2) chain(i,3)]);
  evalChain(i,2) = twoModeSurrogate([chain(i,4) chain(i,5)]);
end
p = 0:0.005:1;
[X Y Z] = twoModeSurrogateGrid(p,p);
pcolor(X,Y,Z);             
shading flat;              
hold on;                   
pointsSizeScale = 60;
scatter(chain(:,2),chain(:,3),pointsSizeScale*weights(:,1),'filled','red');
scatter(chain(:,4),chain(:,5),pointsSizeScale*weights(:,2),'filled','black');
c=colorbar; c.Label.String = 'Y(theta_1,theta_2)';
xlabel('theta_1','FontName',fname,'fontsize',fsize);
ylabel('theta_2','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
hold off;

figure()
surf(X,Y,Z);             
shading flat; 
%alpha 0.5;             
hold on;
scatter3(chain(:,2),chain(:,3),evalChain(:,1),pointsSizeScale*weights(:,1),'filled','red');
scatter3(chain(:,4),chain(:,5),evalChain(:,2),pointsSizeScale*weights(:,2),'filled','black');
xlabel('theta_1','FontName',fname,'fontsize',fsize);
ylabel('theta_2','FontName',fname,'fontsize',fsize);
zlabel('Y(theta_1,theta_2)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
daspect([1 1 4]);
hold off;

figure()
chainElementsIndices = 1:length(chain);
chainElementsIndices = chainElementsIndices + burning;
plot(chainElementsIndices,weights(:,1),'red');
hold on;
plot(chainElementsIndices,weights(:,2),'black');
xlabel('iteration','FontName',fname,'fontsize',fsize);
ylabel('mixture weight','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
hold off;
