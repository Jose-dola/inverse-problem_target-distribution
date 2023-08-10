function err = E(bound1,bound2,theta,numberOfGaussians,parametersDimension,numberOfSamplesTraining ...
                            ,numberOfSamplesTrueDensity,trueDensitySamples,MMDtrueDensitySamples,kernelSigmaSquare)

  [weights,means,sigmas] = mcmcParamToRealParam(theta,numberOfGaussians,parametersDimension);
  
%  samples = twoModeSurrogateMultipleEvaluation ...
%              ( getSamplesFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension,numberOfSamplesTraining) );
%
%
%  err     = MMDtrueDensitySamples ...
%            + ( sum( exp( (-(samples.' - samples).^2) / (2*kernelSigmaSquare) ) , 'all') - numberOfSamplesTraining ) ...
%                                                         / ( numberOfSamplesTraining*(numberOfSamplesTraining-1) )   ...
%            - ( 2*sum( exp( (-(samples.' - trueDensitySamples).^2) / (2*kernelSigmaSquare) ) , 'all') ) ...
%                                                         / ( numberOfSamplesTraining*numberOfSamplesTrueDensity ) ;
%  err = abs(err);
%

  samplesInputSpace = getSamplesFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension,numberOfSamplesTraining);
  samplesSurrogate1 = twoModeSurrogateMultipleEvaluation(samplesInputSpace);
  samplesSurrogate2 = twoModeSurrogate2MultipleEvaluation(samplesInputSpace);

  err = sum(samplesSurrogate2 > bound2) + sum(samplesSurrogate1 > bound1);

end

