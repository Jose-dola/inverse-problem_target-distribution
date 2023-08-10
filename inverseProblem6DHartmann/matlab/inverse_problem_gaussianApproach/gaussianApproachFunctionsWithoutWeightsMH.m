function p = gaussianApproachFunctionsWithoutWeightsMH(theta,ntheta,x)

m = theta( 1          : ntheta/2 );
l = theta( ntheta/2+1 : ntheta   );

aux = x.'-m;
p   = sum((exp(-(aux.*aux)./(2*l))).');

end
