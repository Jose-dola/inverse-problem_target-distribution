function fillHisto2(theta,numberOfGaussians,parametersDimension,h,numberOfSamples)
  h.resetObservations();
  [weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension);
  for i = 1:numberOfSamples
    h.pointToHisto( twoModeSurrogate2( getRandomFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension) ) );
  end
  h.normalize(); 
end

