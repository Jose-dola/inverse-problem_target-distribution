function modelout = modelpolyfun(theta,xres)

% for ii=1:length(xres)
% modelout(ii) = polyapp(xres(ii))*theta.';   % MODEL YOUR FUNCTION
% end

modeloutun = theta*polyapp(xres);
modelout = modeloutun./sum(modeloutun);

end


function p = polyapp(x)

nord = 4;
nx = length(x);
p = zeros(nord,nx);

p(1,:) = 1;
p(2,:) = x;
p(3,:) = (3.*x.^2 - 1)./2; 
p(4,:) = (35.*x.^4 - 30.*x.^2 + 3)./8;

% p=p+1;

end
