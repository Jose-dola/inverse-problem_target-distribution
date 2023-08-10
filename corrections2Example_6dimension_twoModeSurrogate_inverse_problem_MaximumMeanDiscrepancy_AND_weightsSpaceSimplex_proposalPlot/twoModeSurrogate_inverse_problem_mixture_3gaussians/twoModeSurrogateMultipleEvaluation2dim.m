function z = twoModeSurrogateMultipleEvaluation2dim(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% X: matrix whose columns are vectors [x,y]
%      x: x coordinate
%      y: y coordinate
%
% OUTPUT:
% z = vector whose elements are evaluations of -> surrogate(x,y) --- according to X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(X);

mode1 = [1/3 2/3];
mode2 = [2/3 1/3];
sigma = 0.05;
v     = mode1 - mode2;
scaleToOne = 1 + exp(-(v*v.' / sigma));

X1 = X - mode1.';
X2 = X - mode2.';
for i = 1:n                   
  x1(i) = X1(:,i).'*X1(:,i);
  x2(i) = X2(:,i).'*X2(:,i);
end

z = 2*( ( -exp(-x1./sigma) - exp(-x2./sigma) ) ./ scaleToOne );

end
