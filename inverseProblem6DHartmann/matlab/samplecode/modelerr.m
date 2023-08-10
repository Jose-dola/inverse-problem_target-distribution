
function errout = modelerr(theta,tpdf,xres)

nres=length(xres);

% for ii=1:10
% t(ii) = theta(ii);
% end

% errout = sum((tpdf(xres) - modelfun(theta)).^2.)/nres;
normtarget = tpdf(xres)./sum(tpdf(xres));
errout = log(sum((normtarget - modelpolyfun(theta,xres)).^2));

end
