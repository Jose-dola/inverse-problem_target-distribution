function z = twoModeSurrogate2MultipleEvaluation(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% X: matrix whose columns are vectors [x1 ... x10]
%
% OUTPUT:
% z = vector whose elements are evaluations of -> surrogate2(x1 ... x10) --- according to X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%z = zeros(1,size(X,2));
%for i = 0:4
%  z = z+twoModeSurrogateMultipleEvaluation2dim( X(2*i+1:2*i+2,:) );
%end

dim = size(X);

numberOfSamples = dim(2);
parametersDimension = dim(1);

for i=1:numberOfSamples
        z(i) = 1;
        for j=1:parametersDimension
                z(i) = z(i)*X(j,i);
        end
%	z(i) = z(i)^(1/5);
end


end
