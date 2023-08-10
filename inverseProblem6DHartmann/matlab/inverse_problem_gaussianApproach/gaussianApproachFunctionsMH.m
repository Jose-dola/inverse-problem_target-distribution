function p = gaussianApproachFunctionsMH(theta,ntheta,x)

alpha = theta( 1            : ntheta/3   );
m     = theta( ntheta/3+1   : 2*ntheta/3 );
l     = theta( 2*ntheta/3+1 : ntheta     );

exp(l);

aux = x.'-m;
p   = (exp(-(aux.*aux)./(2*l))*(alpha.')).';

end
