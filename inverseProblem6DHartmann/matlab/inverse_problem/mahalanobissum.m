function msum = mahalanobissum(msum_prev,x,mu,sigma)
  msum = msum_prev + ( (x-mu)*(x-mu) / 2*sigma*sigma );
end

