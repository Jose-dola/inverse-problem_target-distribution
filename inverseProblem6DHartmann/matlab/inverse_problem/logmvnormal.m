function l = logmvnormal(x,mu,covmatrix)
  dif = x-mu;
  l = dif*covmatrix*dif.';
end
