function q = proposal(p,sigma)
  s = size(p);
  q = zeros(s);
  for i = 1:s(2)
    aux = normrnd(p(i),sigma);
    while aux < 0 | aux > 1
      aux = normrnd(p(i),sigma);
    end
    q(i) = aux;
  end
end
