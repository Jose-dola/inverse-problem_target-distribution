function x = gradientDescent(x0, f, fgradient, stepSize, MAX_ITER, TOL)
  i = 0;
  x = x0;
  while i < MAX_ITER && f(x) > TOL
    i     = i + 1;
    x     = x - stepSize*fgradient(x);
  end
end
