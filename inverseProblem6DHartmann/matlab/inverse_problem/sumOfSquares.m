function s = sumOfSquares(x, data)
  s = (data.sumOfQuadraticTerms + data.nTerms*x*x - 2*x*data.sumOfTerms) / data.sigma_model;
end
