clear all; close all; clc;
fsize=28; fname='times';

addpath('../posteriorPredictiveApproach','-end')
addpath('../target_pdf','-end')
addpath('../histoStructure','-end')

sigmoid = @(x) exp(x) / (exp(x)+1);

matrixTransform3Dto2D = [-1 1 0; 0 0 1];

dim        = 2;
ngaussians = 3;
m     = [0.5 0.5];
covar = 1e-1*eye(ngaussians-1);

nrandom = 1e5;
Random  = mvnrnd(m,covar,nrandom);

theta = [sigmoid(m(1)) , sigmoid(m(2)) ,1,1,1,1,1,1];
mInTriangle = mcmcParamToRealParam(theta,ngaussians,dim);
mInTriangle = matrixTransform3Dto2D*mInTriangle(1:3).';
for i=1:nrandom
  theta(1)  = sigmoid(Random(i,1));%exp(Random(i,1)) / (exp(Random(i,1)) + 1);
  theta(2)  = sigmoid(Random(i,2));%exp(Random(i,2)) / (exp(Random(i,2)) + 1);
  realParam = mcmcParamToRealParam(theta,ngaussians,dim);
  p(i,:)    = matrixTransform3Dto2D*realParam(1:3).';
end

figure();
plot([-1 0 1].', [0 1 0].');
hold on;
scatter(p(:,1),p(:,2),'.');
scatter(mInTriangle(1),mInTriangle(2)); 
hold off;

npoints = 1e1;
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
    realParam = mcmcParamToRealParam(theta,ngaussians,dim);
    p0        = matrixTransform3Dto2D*realParam(1:3).';
    X(i,j)    = p0(1);
    Y(i,j)    = p0(2);
    n=0;
   % aux = p-p0
    for k=1:nrandom
      aux = p0.'-p(k,:);
      if sqrt(aux*aux.') < r
        n = n+1;
      end
    end
    C(i,j) = n;
  end
end

pcolor(X,Y,C);




