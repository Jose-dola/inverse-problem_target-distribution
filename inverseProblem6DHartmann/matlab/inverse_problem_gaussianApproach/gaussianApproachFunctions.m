function p = gaussianApproachFunctions(x,m,l)

aux = x.'-m;
p   = exp(-0.5*kron(aux.*aux,1./l));


end
