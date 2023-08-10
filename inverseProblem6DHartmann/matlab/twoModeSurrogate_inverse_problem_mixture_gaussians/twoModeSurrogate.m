function z = twoModeSurrogate(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% X: X=[x y]
%      x: x coordinate
%      y: y coordinate
%
% OUTPUT:
% z:: evaluation surrogate(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode1 = [1/3 2/3];
mode2 = [2/3 1/3];
sigma = 0.05;
v     = mode1 - mode2;
scaleToOne = 1 + exp(-(v*v.' / sigma));

X1 = (X(1)-mode1(1)).^2;
Y1 = (X(2)-mode1(2)).^2;
X2 = (X(1)-mode2(1)).^2;
Y2 = (X(2)-mode2(2)).^2;

z = 2*( ( -exp((-X1-Y1)/sigma) - exp((-X2-Y2)/sigma) ) / scaleToOne ); 

end
