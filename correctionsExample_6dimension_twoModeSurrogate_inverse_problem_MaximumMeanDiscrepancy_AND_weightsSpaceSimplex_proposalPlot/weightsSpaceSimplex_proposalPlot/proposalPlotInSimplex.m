clear all; close all; clc;
fsize=28; fname='times';

addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
addpath('../histoStructure','-end')

sigmoid    = @(x) exp(x) / (exp(x)+1);
invSigmoid = @(x) log(x/(1-x));

matrixTransform3Dto2D = [-1 1 0; 0 0 1];

dim        = 2;
ngaussians = 3;
%m     = [0 log(2)];
m     = [invSigmoid(0.5) invSigmoid(0.1)];
covar = 1e-1*eye(ngaussians-1);

nrandom = 1e4;
Random  = mvnrnd(m,covar,nrandom);

theta = [sigmoid(m(1)) , sigmoid(m(2))];
mInTriangle = weightsParamToWeightsInSimplex(theta,ngaussians,dim);
mInTriangle = matrixTransform3Dto2D*mInTriangle(1:3).';
for i=1:nrandom
  theta(1)  = sigmoid(Random(i,1));
  theta(2)  = sigmoid(Random(i,2));
  realParam = weightsParamToWeightsInSimplex(theta,ngaussians,dim);
  p(i,:)    = matrixTransform3Dto2D*realParam(1:3).';
end

figure();
plot([-1 0 1].', [0 1 0].');
hold on;
scatter(p(:,1),p(:,2),'.');
scatter(mInTriangle(1),mInTriangle(2),'*','black'); 

npoints = 1e2;
x = linspace(0,1,npoints);
%for i=1:npoints
%  X(:,i) = x;
%end
%Y = X';

r = 1e-2;
for i=1:npoints
  for j=1:npoints
    theta(1)  = x(i);
    theta(2)  = x(j);
    realParam = weightsParamToWeightsInSimplex(theta,ngaussians,dim);
    p0        = matrixTransform3Dto2D*realParam(1:3).';
    X(i,j)    = p0(1);
    Y(i,j)    = p0(2);
    C(i,j) = sum( sum((p-p0.').^2,2) < r );
%    n=0;
%    for k=1:nrandom
%      aux = p0.'-p(k,:);
%      if sqrt(aux*aux.') < r
%        n = n+1;
%      end
%    end
%    C(i,j) = n;
  end
end

figure()
pcolor(X,Y,C);
shading flat; 



