function m = mahalanobisDistance(x,sample,sigma)
  aux = arrayfun(@(y) mahalanobis(y,x,sigma), sample);
  m   = sum(aux);
  %s = size(sample);
  %nTerms = s(2);
  %sumOfTerms = sum(sample);
  %sumOfQuadraticTerms = sum(sample.*sample);
  %mm   = mahalanobisDistanceToSample(x, nTerms, sumOfTerms, sumOfQuadraticTerms, sigma);
  %m-mm < 1e-14 & -1e-14 < m-mm
end

