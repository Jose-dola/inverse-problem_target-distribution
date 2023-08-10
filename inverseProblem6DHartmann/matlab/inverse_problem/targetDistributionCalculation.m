%%% Finding displaced log-normal pdf with maximum in 'Xmax' 
%%% and a tail such that P(x>'Xtail')='tailWeight'.
%%% The evaluation in values less than 'Xzero' is 0.
%%% Then the code generates 'nSample' samples of this distribution.
%%%

global Xmax;
global tailWeight;
global Xzero;
global Xtail;
nSample    = 1e2;
Xmax       = -3.322;
tailWeight = 0.1;
Xzero      = Xmax - 0.1;
Xtail      = -2;

x = [0; 0.1];
while norm(f(x)) > 1e-12
    x = x - Jacobian(x) \ (f(x).');
end

mu    = x(1);
sigma = x(2);

xPlot = Xzero+1e-14:1e-2:1;
%yPlot = arrayfun(@(x) logNormal(x,mu,sigma,Xzero), xPlot);
yPlot = lognpdf(xPlot-Xzero,mu,sigma);
%
sample = random('Lognormal',mu,sigma,[1 nSample]) + Xzero;
%edges  = [Xzero-0.2:1e-3:2]
%hist   = histogram(sample,edges);
%
plot(xPlot,yPlot);
%hist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% posterior predictive
lGrid     = 1e-1;
%ranges    = cell(6,1);
%for i = 1:6
%    ranges{i} = 0:lGrid:1;
%end
%outputs = cell(size(ranges));
%display([ranges{:}]);
%[outputs{:}] = ndgrid(ranges{:});
%outputs{:}

aux = 0:lGrid:1;
s   = size(aux);
n   = s(2)^6;
parametersGrid = zeros([n 6]);
i = 1;
for x1 = aux
    for x2 = aux
        for x3 = aux
            for x4 = aux
                for x5 = aux
                    for x6 = aux
                        parametersGrid(i,:)=[x1 x2 x3 x4 x5 x6];
                        i = i + 1;
                    end
                end
            end
        end
    end
end

hartmanPoints = zeros([1 n]);
for i = 1:n
    hartmanPoints(i) = hart6( parametersGrid(i,:) );
end

edges  = [-3:1e-3:0];
figure()
histogram(hartmanPoints,edges);

%parameter  = [0.4 0.2 0.2 0.4 0.2 0.2];
%h          = hart6(parameter);
sigma      = 1e-2;

hartmanProbs = zeros([1 n]);
for i = 1:n
  aux = arrayfun(@(y) mahalanobis(y,hartmanPoints(i),sigma), sample);
  hartmanProbs(i) = sum(aux);
end

figure()
plot(hartmanPoints, exp(-hartmanProbs)./sum(exp(-hartmanProbs)), '.b'); hold on;
histogram(hartmanPoints,edges,'Normalization','pdf');

