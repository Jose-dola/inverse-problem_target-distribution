N = 6;
a = rand(1,N).*rand(1,N).*rand(1,N);
t = 0.01;
f = @(x) gaussians(x,a,N,t);
x = linspace(-5,5,1e3);
fx = arrayfun(f,x);
plot(x,fx);

function f = gaussians(x,a,N,t)
  f = 0;
  for n = 1:N
    aux = x-n*t;
    f = f + a(n)*exp(-(aux*aux));
  end
end
