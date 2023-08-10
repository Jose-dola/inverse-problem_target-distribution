function m = mahalanobis(x,mu,sigma)
  m = (x-mu)*(x-mu) / 2*sigma*sigma;
end
