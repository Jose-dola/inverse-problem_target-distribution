%% Data
randn('seed',1)                            % Used for reproducibility
data = [randn(100,1)-10; randn(100,1)+10]; % Two Normals mixed

% True
phi  = @(x) exp(-.5*x.^2)/sqrt(2*pi);      % Normal Density
tpdf = @(x) phi(x+10)/2+phi(x-10)/2;       % True Density

% Kernel
h      = std(data)*(4/3/numel(data))^(1/5); % Bandwidth estimated by Silverman's Rule of Thumb
kernel = @(x) mean(phi((x-data)/h)/h);      % Kernel Density
kpdf   = @(x) arrayfun(kernel,x);           % Elementwise application

% Plot
figure(1), clf, hold on
x = linspace(-25,+25,1000);                % Linear Space
plot(x,tpdf(x))                            % Plot True Density
plot(x,kpdf(x))                            % Plot Kernel Density with Silverman's Rule of Thumb
% kde(data)                                  % Plot Kernel Density with Solve-the-Equation Bandwidth


% Plot
figure(2), clf, hold on
x = linspace(-25,+25,1000);                % Linear Space
semilogy(x,tpdf(x))                            % Plot True Density
semilogy(x,kpdf(x))                            % Plot Kernel Density with Silverman's Rule of Thumb



%% experiment

clear all; close all; clc;
addpath('../mcmcstat','-end')
fsize=18; fname='times';

nres  = 100;
xres  = linspace(-1,1,nres);
nord  = 4;
sigma = 0.3.^2;

%phi  = @(x) exp(-.5*x.^2./2/sigma)/sqrt(2*pi*sigma); % Normal Density
phi  = @(x) exp(-x.^2./2/sigma)/sqrt(2*pi*sigma); % Normal Density
tpdf = @(x) (phi(x-0.5)+phi(x+0.5))./2;           % Target Density
% polyapp = @(x) legendreP([1:nord],x);

figure()
tpdfOVERxres = tpdf(xres);
plot(xres,tpdfOVERxres./sum(tpdfOVERxres))

% modelfun = @(theta) arrayfun(polyapp,xres)*theta(:).';   % MODEL YOUR FUNCTION
% ssfun = @(theta,dummy) modelerr(theta,modelfun);
ssfun = @(theta,dummy) modelerr(theta,tpdf,xres);

for jj = 1:5
  tuni = -1+2.*rand(1,nord);
  -log(sum((tpdfOVERxres - modelpolyfun(tuni,xres)).^2)./nres)
end



params = {
    {'theta1', 0.5, -10, 10}
    {'theta2', 0.5, -10, 10}
    {'theta3', 0.5, -10, 10}
    {'theta4', 0.5, -10, 10}
%     {'theta5', 0.5, -10, 10}
%     {'theta6', 0.5, -10, 10}
%     {'theta7', 0.5, -10, 10}
%     {'theta8', 0.5, -10, 10}
%     {'theta9', 0.5, -10, 10}
%     {'theta10', 0.5, -10, 10}
    };

model.ssfun  = ssfun;
model.sigma2 = 0.01^2;
options.method = 'am';
options.nsimu = 100000;
options.updatesigma = 0;
options.qcov = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
% options.adaptint = 10000;

[res0,chain0,s2chain0] = mcmcrun(model,[],params,options);

burnin = 1000; thin = 10;
chain = chain0(burnin:thin:end,:); s2chain = s2chain0(burnin:thin:end,:);

figure(2); clf
mcmcplot(chain,[],res0,'chainpanel');
figure(3); clf
mcmcplot(chain,[],res0,'pairs');

for ii=1:size(chain,1)
    errchain(ii) = ssfun(chain(ii,:),[]);
end
postpred = modelpolyfun(mean(chain,1),xres);
figure()
plot(errchain);

figure()
plot(xres,postpred, '-r'); hold on;
plot(xres,tpdf(xres)./sum(tpdf(xres)), '-b')
% plot(1:nres,modelpolyfun(0.5.*ones(1,nord),xres), '-k')
hold off


% figure()
% plot(chain(:,1),chain(:,2),'-b')

%%


%n=1e6;
%
%x0 = randn(n,1);
%
%for ii=1:n
%
%    xm(ii) = randn(1);
%    
%    while (xm(ii) < -0.8) || (xm(ii) > 0.8)
%        xm(ii) = randn(1);
%    end    
%end
%
%figure()
%histogram(x0, 'Normalization', 'pdf'); hold on;
%histogram(xm, 'Normalization', 'pdf'); hold off;
%grid on






