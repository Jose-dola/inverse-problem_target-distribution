clear all;

n = 5000;
for i = 1:n
  p    = randomPoint1();
  x(i) = p(1);
  y(i) = p(2);
end
figure()
plot(x,y,'.');
hold off;

for i = 1:n
  p    = randomPoint2();
  x(i) = p(1);
  y(i) = p(2);
end
figure()
plot(x,y,'.');
hold off;

for i = 1:n
  p    = randomPoint2WithSQRT();
  x(i) = p(1);
  y(i) = p(2);
end
figure()
plot(x,y,'.');
hold off;

for i = 1:n
  p    = randomPointDirichlet();
  x(i) = p(1);
  y(i) = p(2);
end
figure()
plot(x,y,'.');
hold off;

function p = randomPoint1()
  mu1 = rand();
  mu2 = mu1+rand()*(1-mu1);
  p   = [0 sqrt(6)/2] - mu1*sqrt(2)*[1 0] + mu2*0.5*sqrt(2)*[1 -sqrt(3)];
end

function p = randomPoint2()
  mu2 = rand();
  mu1 = mu2*rand();
  p   = [0 sqrt(6)/2] - mu1*sqrt(2)*[1 0] + mu2*0.5*sqrt(2)*[1 -sqrt(3)];
end

function p = randomPoint2WithSQRT()
  mu2 = sqrt( rand() );
  mu1 = mu2*rand();
  p   = [0 sqrt(6)/2] - mu1*sqrt(2)*[1 0] + mu2*0.5*sqrt(2)*[1 -sqrt(3)];
end

function p = randomPointDirichlet()
  A = [-sqrt(2)/2 sqrt(2)/2 0; 0 0 sqrt(6)/2];
  for i = 1:3
    x(i) = gamrnd(1,1);
  end
  x = x./sum(x);
  p = A*x.';
end
