function [weights,means,sigmas] = weightsParamToWeightsInSimplex(weightsParam,numberOfGaussians,parametersDimension)
  aux  = weightsParam(numberOfGaussians-1);
  weights(numberOfGaussians) = 1-aux;
  for i = 1:numberOfGaussians-2
    aux2       = weightsParam(numberOfGaussians-i-1)*aux;
    weights(numberOfGaussians-i) = aux - aux2;
    aux        = aux2; 
  end
  weights(1) = aux;
end

