function x = newton( f, x0, tol, nmax )
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % NEWTON'S METHOD TO FIND THE ROOT OF f
  % JACOBIAN OF f IS CALCULATED USING A NUMERICAL APPROXIMATION
  % param:
  %   f    : function
  %   x0   : initial point
  %   tol  : tolerance ( if norm(f(x)) < tol --> return x )
  %   nmax : maximum number of iterations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n  = length(x0);
  k  = 1;
  x  = x0
  fx = f(x); 
  while (norm(fx) >= tol) && (k <= nmax)
    x  = x - ( Jacobian(f,x,n) \ (fx.') ).'
    fx = f(x);
    k  = k+1;
  end
end

function df = Jacobian(f,x,n)
  delta = 1e-8*norm(x);
  for i = 1:n
    xgreater = x;
    xlesser  = x;
    xgreater(i) = xgreater(i)+delta;
    xlesser(i)  = xlesser(i)-delta;
    df(:,i) = (f(xgreater)-f(xlesser))./(2*delta);
  end
end

