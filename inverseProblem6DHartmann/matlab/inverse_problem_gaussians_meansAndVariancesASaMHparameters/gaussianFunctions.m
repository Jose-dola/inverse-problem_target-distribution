function p = gaussianFunctions(x,m,l)

aux = x.'-m;
p   = exp(-(aux.*aux)./(2*l));

end
