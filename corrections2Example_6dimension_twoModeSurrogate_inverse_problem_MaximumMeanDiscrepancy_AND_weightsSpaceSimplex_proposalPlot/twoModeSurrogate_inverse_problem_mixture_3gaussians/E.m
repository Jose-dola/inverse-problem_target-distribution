function err = E(theta,numberOfGaussians,parametersDimension,numberOfSamplesTraining ...
                            ,truePDFmiddlePoints,histoStructure,kernelSigmaSquare)

  [weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension);
  
  samples = twoModeSurrogateMultipleEvaluation ...
              ( getSamplesFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension,numberOfSamplesTraining) );

  histoStructure.resetObservations();
  for i=1:length(samples)
	histoStructure.pointToHisto(samples(i));
  end
  histoStructure.normalize();
  err = histoStructure.L2error(truePDFmiddlePoints);  

end

