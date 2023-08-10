function J = Jacobian(x)
  mu       = x(1);
  sigma    = x(2);
  delta = 1e-8*norm([mu,sigma]);

  dfdmu    = (f([mu+delta,sigma])-f([mu-delta,sigma]))/(2*delta);
  dfdsigma = (f([mu,sigma+delta])-f([mu,sigma-delta]))/(2*delta);

  J = [dfdmu; dfdsigma].';
end
