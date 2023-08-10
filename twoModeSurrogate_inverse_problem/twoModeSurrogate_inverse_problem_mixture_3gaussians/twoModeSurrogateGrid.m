function [X Y Z] = twoModeSurrogateGrid(x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% x: x coordinates
% y: y coordinates
%
% OUTPUT:
% X: x coordinate grid
% Y: y coordinate grid
% Z: evaluation surrogate(x,y) over the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X Y] = meshgrid(x,y);

mode1 = [1/3 2/3];
mode2 = [2/3 1/3];
sigma = 0.05;
v     = mode1 - mode2;
scaleToOne = 1 + exp(-(v*v.' / sigma));

X1 = (X-mode1(1)).^2;
Y1 = (Y-mode1(2)).^2;
X2 = (X-mode2(1)).^2;
Y2 = (Y-mode2(2)).^2;

Z = 2*( ( -exp((-X1-Y1)/sigma) - exp((-X2-Y2)/sigma) ) / scaleToOne ); 

end
