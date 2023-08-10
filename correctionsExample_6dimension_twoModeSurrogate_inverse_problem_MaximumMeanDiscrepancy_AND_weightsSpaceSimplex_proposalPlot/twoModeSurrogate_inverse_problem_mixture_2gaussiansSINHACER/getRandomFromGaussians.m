function r = getRandomFromGaussians(weights,means,sigmas,numberOfGaussians,parametersDimension)
  rd = rand();
  w  = 0;
  i  = 0;
  while w < rd
    w = w + weights(i+1);
    i = i+1;
  end
  r = mvnrnd(means(i,:),sigmas(i)*eye(parametersDimension));
end

