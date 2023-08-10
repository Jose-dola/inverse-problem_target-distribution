function [weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension)
  weightsmcmc = [theta(1:numberOfGaussians-1) 1];
  weights     = zeros(1,numberOfGaussians);
  means       = zeros(numberOfGaussians,parametersDimension);
  means(1,:) = theta( numberOfGaussians : numberOfGaussians+parametersDimension-1 );
  firstIndex = numberOfGaussians+parametersDimension;
  for i = 2:numberOfGaussians
    means(i,:) = theta( firstIndex : firstIndex+parametersDimension-1 );
    firstIndex = firstIndex+parametersDimension;
    j = numberOfGaussians-i+1;
    weightsmcmc(j) = ( weightsmcmc(j)*weightsmcmc(j+1) )^(1/j);
    weights(j+1)   = weightsmcmc(j+1)-weightsmcmc(j);
  end
  weights(1) = weightsmcmc(1);
  sigmas = theta(firstIndex : end);
end

