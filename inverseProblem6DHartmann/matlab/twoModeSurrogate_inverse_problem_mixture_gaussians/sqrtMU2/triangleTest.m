clear all;

n = 5000;

for i = 1:n
  p    = randomPoint();
  x(i) = p(1);
  y(i) = p(2);
end
figure()
plot(x,y,'.');
hold off;

function p = randomPoint()
  A = [-sqrt(2)/2 sqrt(2)/2 0; 0 0 sqrt(6)/2];
  x = [rand(1,2) 1 1 1 1 1 1];
  [weights,means,sigmas] = mcmcParamToRealParam(x,3,1);
  p = A*weights.';
end
