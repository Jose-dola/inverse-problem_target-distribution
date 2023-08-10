function z = twoModeSurrogateMultipleEvaluation(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% X: matrix whose columns are vectors [x1 ... x10]
%
% OUTPUT:
% z = vector whose elements are evaluations of -> surrogate(x1 ... x10) --- according to X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = zeros(1,size(X,2));
for i = 0:4
  z = z+twoModeSurrogateMultipleEvaluation2dim( X(2*i+1:2*i+2,:) );
end

end
