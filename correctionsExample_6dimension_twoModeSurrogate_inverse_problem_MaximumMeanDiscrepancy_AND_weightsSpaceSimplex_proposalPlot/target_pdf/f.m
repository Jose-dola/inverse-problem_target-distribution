function f = f(x,Xmax,Xzero,Xtail,tailWeight)
  f     = [0 0];
  mu    = x(1);
  sigma = x(2);
  
  f(1) = mvncdf(log(Xtail-Xzero),mu,sigma) - 1 + tailWeight;
  f(2) = log(Xmax-Xzero) + sigma*sigma - mu;
end
