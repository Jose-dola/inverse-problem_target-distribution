function p = gaussianApproachFunctionsSingleVarianceMH(theta,nalpha,x)

t     = theta(1);
alpha = theta(2:end);
for i = 1:nalpha
  m(i) = (i-1)*t;
end

aux = x.'-m;
p   = (exp(-(aux.*aux))*(alpha.')).';

end
