function J = Jacobian(x,Xmax,Xzero,Xtail,tailWeight)
  mu       = x(1);
  sigma    = x(2);
  delta    = 1e-8*norm([mu,sigma]);

  dfdmu    = (f([mu+delta,sigma],Xmax,Xzero,Xtail,tailWeight)-f([mu-delta,sigma],Xmax,Xzero,Xtail,tailWeight))/(2*delta);
  dfdsigma = (f([mu,sigma+delta],Xmax,Xzero,Xtail,tailWeight)-f([mu,sigma-delta],Xmax,Xzero,Xtail,tailWeight))/(2*delta);

  J = [dfdmu; dfdsigma].';
end
