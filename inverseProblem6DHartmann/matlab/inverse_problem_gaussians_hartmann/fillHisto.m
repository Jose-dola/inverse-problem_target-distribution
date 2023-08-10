function fillHisto(theta,numberOfGaussians,parametersDimension,h,numberOfSamples)
  h.resetObservations();
  [weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension);
  for i = 1:numberOfSamples
    h.pointToHisto( hart6( getRandomFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension) ) );
  end
  h.normalize(); 
end

