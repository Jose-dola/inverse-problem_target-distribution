clear all; close all; clc;
fsize=28; fname='times';

addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
addpath('../histoStructure','-end')
% Number of gaussians
ngaussians = 5;
% Parameters Dimension
parametersDimension = 10;
%Gaussian kernel sigma
kernelSigmaSquare = 0.001;
%% Target Density
Xmax       = -1.9*5-0.2;
Xzero      = Xmax - 0.2;
Xtail      = -7.5;
tailWeight = 0.2;
%[tpdf, randtpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight);

Xzero      = -10;
truepdfsigma = 1;
truepdfmean  = -5.5;
randtpdf  = @(x) randn(x)*truepdfsigma + truepdfmean;
tpdf      = @(x) normpdf(x,truepdfmean,truepdfsigma);

%figure()
%histogram(randtpdf([1,5000]),'Normalization','probability'); hold on;
%x = linspace(-10,2,1000);
%plot(x,tpdf(x))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histo Structure
h = histo(Xzero,1,200);
%% true density
trueDensity = tpdf( h.intervalrepresentativePoints() );
%numberOfSamplesTrueDensity = 1e2;
%trueDensitySamples = randtpdf([1 numberOfSamplesTrueDensity]);
%MMDtrueDensitySamples = ( sum( exp( (-(trueDensitySamples.' - trueDensitySamples).^2) / (2*kernelSigmaSquare) ) , 'all') - numberOfSamplesTrueDensity ) ...
%                                                         / ( numberOfSamplesTrueDensity*(numberOfSamplesTrueDensity-1) );

%% mcmc model error
numberOfSamplesTraining = 5000;
%boundSurrogate1 = -4;
%boundSurrogate2 = 0;
modelerr = @(theta) E( ...
	           theta,ngaussians,parametersDimension,numberOfSamplesTraining ...
                                         ,trueDensity,h,kernelSigmaSquare);
modelSigma2 = 0.0003;
modelPriorSigma2S20 = modelSigma2;
modelPriorSigma2N0  = 1e2;
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
options.nsimu       = 1e6;%1e6;
options.updatesigma = 1;
options.qcov        = proposalCov;
burning             = 1;%10000;
thin                = 1;
%% POSTERIOR PREDICTIVE APPROACH -- RUNING MCMC
[res0, chain, s2chain, sschain] = posteriorPredictiveApproach(modelerr, modelSigma2, modelPriorSigma2S20, modelPriorSigma2N0, params, options, burning, thin, false);

%% RESULTS %%
% saved the chain
fi = fopen('chain.doubleMatrix.bin','w');
fwrite(fi, chain, 'double');   
fclose(fi); 
fi = fopen('chainSize.doubleVector.bin','w');
fwrite(fi, size(chain), 'double');   
fclose(fi); 

% errors and best approximation in the chain
errchain = sschain.';
j = 1;
for i=2:size(chain,1)
  if (errchain(i) < errchain(j)) j = i; end
end
bestChain = chain(j,:);
fi = fopen('bestChain.double.bin','w');
fwrite(fi, bestChain, 'double'); 
fclose(fi);  

% plot chain error
f = figure();
plot(errchain);
savefig(f,'errchain.fig');
%print(f,'errchain','-depsc');


% plot true density and predicted densities for different theta values
f1 = figure();
lineWidth = 2;

meanChain = mean(chain,1);
fi = fopen('meanChain.double.bin','w');
fwrite(fi, meanChain, 'double'); 
fclose(fi);  
%5fillHisto(meanChain,ngaussians,parametersDimension,h,1e3);
%[f,xi] = h.ksdens;
%plot(xi,f,'-r','LineWidth',lineWidth);

hold on;

lastChain = chain(end,:);
fi = fopen('lastChain.double.bin','w');
fwrite(fi, lastChain, 'double'); 
fclose(fi);  
%fillHisto(lastChain,ngaussians,parametersDimension,h,1e3);
%[f,xi] = h.ksdens;
%plot(xi,f,'-g','LineWidth',lineWidth);

fillHisto(theta0,ngaussians,parametersDimension,h,1e3);
[f,xi] = h.ksdens;
plot(xi,f,'-r','LineWidth',lineWidth);

fillHisto(bestChain,ngaussians,parametersDimension,h,1e3);
[f,xi] = h.ksdens;
plot(xi,f,'-k','LineWidth',lineWidth);

plot(h.intervalrepresentativePoints,trueDensity, '-b','LineWidth',lineWidth);

xlabel('y','FontName',fname,'fontsize',fsize);
ylabel('P(y)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');

hold off;
savefig(f1,'surrogate1histo.fig');
%print(f1,'surrogate1histo','-depsc');

%f2 = figure();
%lineWidth = 2;
%
%fillHisto2(meanChain,ngaussians,parametersDimension,h,1e3);
%[f,xi] = h.ksdens;
%plot(xi,f,'-r','LineWidth',lineWidth);
%
%hold on;
%
%fillHisto2(lastChain,ngaussians,parametersDimension,h,1e3);
%[f,xi] = h.ksdens;
%plot(xi,f,'-g','LineWidth',lineWidth);
%
%fillHisto2(theta0,ngaussians,parametersDimension,h,1e3);
%[f,xi] = h.ksdens;
%plot(xi,f,'-c','LineWidth',lineWidth);
%
%fillHisto2(bestChain,ngaussians,parametersDimension,h,1e3);
%[f,xi] = h.ksdens;
%plot(xi,f,'-k','LineWidth',lineWidth);
%
%xline(boundSurrogate2);
%%plot(h.intervalrepresentativePoints,trueDensity, '-b','LineWidth',lineWidth);
%
%xlabel('y','FontName',fname,'fontsize',fsize);
%ylabel('P2(y)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
%
%hold off;
%savefig(f2,'surrogate2histo.fig');
%print(f2,'surrogate2histo','-depsc');

% plot how the means move
color1 = [0.4660 0.6740 0.1880]; % dark green
color2 = [0.6350 0.0780 0.1840]; % dark red
color3 = [0.4940 0.1840 0.5560]; % dark magenta
%
%figure()
for i = 1:length(chain)
  [weights(i,:), means(i,:,:),sigmas(i,:)]   = mcmcParamToRealParam(chain(i,:),ngaussians,parametersDimension);
%  evalChain(i,1) = twoModeSurrogate([chain(i,3) chain(i,4)]);
%  evalChain(i,2) = twoModeSurrogate([chain(i,5) chain(i,6)]);
%  evalChain(i,3) = twoModeSurrogate([chain(i,7) chain(i,8)]);
end
%p = 0:0.005:1;
%[X Y Z] = twoModeSurrogateGrid(p,p);
%pcolor(X,Y,Z);             
%shading flat;              
%hold on;                   
%pointsSizeScale = 60;
%scatter(chain(:,3),chain(:,4),pointsSizeScale*weights(:,1),color1,'filled'); 
%scatter(chain(:,5),chain(:,6),pointsSizeScale*weights(:,2),color2,'filled'); 
%scatter(chain(:,7),chain(:,8),pointsSizeScale*weights(:,3),color3,'filled'); 
%c=colorbar; c.Label.String = 'Y(theta_1,theta_2)';
%xlabel('theta_1','FontName',fname,'fontsize',fsize);
%ylabel('theta_2','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
%hold off;
%
%figure()
%surf(X,Y,Z);             
%shading flat; 
%%alpha 0.5;             
%hold on;
%scatter3(chain(:,3),chain(:,4),evalChain(:,1),pointsSizeScale*weights(:,1),color1,'filled'); 
%scatter3(chain(:,5),chain(:,6),evalChain(:,2),pointsSizeScale*weights(:,2),color2,'filled'); 
%scatter3(chain(:,7),chain(:,8),evalChain(:,3),pointsSizeScale*weights(:,3),color3,'filled'); 
%xlabel('theta_1','FontName',fname,'fontsize',fsize);
%ylabel('theta_2','FontName',fname,'fontsize',fsize);
%zlabel('Y(theta_1,theta_2)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
%daspect([1 1 4]);
%alpha 0.4
%scatter3([bestChain(3) bestChain(5) bestChain(7)].', [bestChain(4) bestChain(6) bestChain(8)].', ...
%        [twoModeSurrogate([bestChain(3) bestChain(4)]) twoModeSurrogate([bestChain(5) bestChain(6)]) twoModeSurrogate([bestChain(7) bestChain(8)])].', ...
%        [500 500 500].','filled','black');
%scatter3([lastChain(3) lastChain(5) lastChain(7)].', [lastChain(4) lastChain(6) lastChain(8)].', ...
%        [twoModeSurrogate([lastChain(3) lastChain(4)]) twoModeSurrogate([lastChain(5) lastChain(6)]) twoModeSurrogate([lastChain(7) lastChain(8)])].', ...
%        [500 500 500].','filled','green');
%hold off;
%
f3 = figure();
chainElementsIndices = 1:length(chain);
chainElementsIndices = chainElementsIndices + burning;
for i=1:ngaussians
  plot(chainElementsIndices,weights(:,i));%,'Color',color1); 
  hold on;
end
%plot(chainElementsIndices,weights(:,2),'Color',color2); 
%plot(chainElementsIndices,weights(:,3),'Color',color3); 
xlabel('iteration','FontName',fname,'fontsize',fsize);
ylabel('mixture weight','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
hold off;
savefig(f3,'weights.fig');
%print(f3,'weights','-depsc');

f4 = figure();
tiledlayout(3,2)
for j=1:ngaussians
  nexttile
  for i=1:parametersDimension
    plot(chainElementsIndices,means(:,j,i)); hold on;
  end
  xlabel('iteration','FontName',fname,'fontsize',fsize);
  ylabel(['mean ',num2str(j)],'FontName',fname,'fontsize',fsize);
  ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
  hold off;
end
savefig(f4,'means.fig');
%print(f4,'means','-depsc');

%close all;
