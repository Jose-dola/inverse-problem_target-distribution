function err = E(theta,numberOfGaussians,parametersDimension,h,numberOfSamples,trueDensity)
  fillHisto(theta,numberOfGaussians,parametersDimension,h,numberOfSamples)
  err = h.L2error(trueDensity);
end

