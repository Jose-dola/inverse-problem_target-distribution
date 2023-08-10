function p = legendreApproachFunctionsRecursive(x,npar)

nx = length(x);
p = zeros(nx,npar);

p(:,1)  = 1;
p(:,2)  = x.';
for n = 3:npar
  p(:,n) = ((2*n-1)/n).*x.*p(:,n-1).' - ((n-1)/n).*p(:,n-2).';
end

end
