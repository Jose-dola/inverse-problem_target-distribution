function y = logNormal(x,mu,sigma,displacement)
  y = exp( -( log(x-displacement) - mu )^2 / (2*(sigma^2)) ) / ( (x-displacement)*sigma*sqrt(2*pi) );
end
