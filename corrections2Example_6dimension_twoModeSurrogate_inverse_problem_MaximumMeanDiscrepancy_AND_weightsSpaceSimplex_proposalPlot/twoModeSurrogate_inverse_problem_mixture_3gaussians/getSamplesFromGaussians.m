function p = getSamplesFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension,numberOfSamples)

for i=1:numberOfSamples
  p(:,i) = getRandomFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension).';
end

%  rd = sort( rand([1 numberOfSamples]) );
%  p  = randn(parametersDimension,numberOfSamples);
%  cumsumweights = cumsum( weights(1:numberOfGaussians-1) );
%  j=1;
%  for i=1:numberOfGaussians-1
%    aux = j;
%    while  j <= numberOfSamples & rd(j) < cumsumweights(i)
%      j = j + 1;
%    end
%    cardinalAccordingToWeight(i) = j - aux;
%    if j>numberOfSamples 
%      cardinalAccordingToWeight(i+1:numberOfGaussians) = 0;
%      break; 
%    end
%  end
%  if j <= numberOfSamples
%    cardinalAccordingToWeight(numberOfGaussians) = numberOfSamples - j + 1;
%  end
%
%  aux = 1;
%  for i=1:numberOfGaussians
%    p(:,aux:aux+cardinalAccordingToWeight(i)-1) = sigmas(i)*p(:,aux:aux+cardinalAccordingToWeight(i)-1) + means(i,:).';
%    aux = aux + cardinalAccordingToWeight(i);
%  end

end

