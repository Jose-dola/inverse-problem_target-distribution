function [pdf, randpdf] = logNormalGenerator(Xmax, Xzero, Xtail, tailWeight)
%%% This function finds a displaced log-normal PDF with the 
%%% maximum in 'Xmax' and the tail such that P(x>'Xtail')='tailWeight'.
%%% The probability in values less than 'Xzero' is 0.
%%%
%%% return: 
%%%   [pdf, randpdf] (vector of two functions):
%%%     pdf     : probability density function of the 
%%%               desired displaced log-normal.
%%%               The only parameter of this function is the 
%%%               point (or points) you want to evaluate.
%%%     randpdf : function that generates random numbers 
%%%               (i.e. samples) acording to the desired 
%%%               displaced log-normal.
%%%               The only parameter of this function is the
%%%               dimension of the sample that is wanted to generate.
%%%               For example, "randpdf([3 1000])" generates a 
%%%               matrix of 3 rows and 1000 columns with random
%%%               numbers generated acording to the desired 
%%%               displaced log-normal.
%%%

aux = [0; 0.1];
while norm(f(aux,Xmax,Xzero,Xtail,tailWeight)) > 1e-12
    aux = aux - Jacobian(aux,Xmax,Xzero,Xtail,tailWeight) \ (f(aux,Xmax,Xzero,Xtail,tailWeight).');
end

mu    = aux(1);
sigma = aux(2);

pdf     = @(x) lognpdf(x-Xzero,mu,sigma);
randpdf = @(x) random('Lognormal',mu,sigma,x) + Xzero;

end
