function m = mahalanobisDistanceToSample(x, nTerms, sumOfTerms, sumOfQuadraticTerms, sigma )
  m = ( sumOfQuadraticTerms + nTerms*x*x - 2*x*sumOfTerms ) / 2*sigma*sigma;
end
