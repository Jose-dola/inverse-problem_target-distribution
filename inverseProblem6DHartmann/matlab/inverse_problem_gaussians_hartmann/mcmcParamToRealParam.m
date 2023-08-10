function [weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension)
  weightsmcmc = [0 theta(1:numberOfGaussians-1) 1];
  weights     = zeros(1,numberOfGaussians);
  means       = zeros(numberOfGaussians,parametersDimension);
  firstIndex  = numberOfGaussians;
  for i = 1:numberOfGaussians
    means(i,:)       = theta( firstIndex : firstIndex+parametersDimension-1 );
    firstIndex       = firstIndex+parametersDimension;
    weightsmcmc(i+1) = weightsmcmc(i)+weightsmcmc(i+1)*(1-weightsmcmc(i));
    weights(i)       = weightsmcmc(i+1)-weightsmcmc(i);
  end
  sigmas = theta(firstIndex : end);
end

