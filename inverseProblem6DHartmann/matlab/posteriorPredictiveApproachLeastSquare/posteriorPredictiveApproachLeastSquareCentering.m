function [postpred,err] = posteriorPredictiveApproachLeastSquareCentering(lowest, highest, n, tpdf, modelApproachFunctions, plots)
%% Points (discretize integral)
x = linspace(lowest,highest,n);
%% Target density evaluations
xcentered = (x-lowest).*(2/(highest-lowest))+lowest;
tpdfEval  = tpdf(x).';
%% model approach functions evaluations
modelApproachFunEval  = modelApproachFunctions(xcentered);
%predictedTheta       = lsqlin(modelApproachFunEval,normalizedTrueDensity,[],[],[],[],lparam,uparam);
%% Least Square Soluction
leastSquareSol = modelApproachFunEval \ tpdfEval;
%% Best model weights that fit the target density
predictedTheta = (leastSquareSol.*n)./( (sum(modelApproachFunEval*leastSquareSol)).*(highest-lowest) );
%% Posterior Predictive PDF
postpred       = @(p) modelApproachFunctions((p-lowest).*(2/(highest-lowest))+lowest)*predictedTheta;
%% Posterior Predictive PDF evaluations
postPredEval   = modelApproachFunEval*predictedTheta;
%% Posterior Predictive PDF error
verror         = tpdfEval-postPredEval;
err            = sqrt( ( (highest-lowest)/n )*( verror.'*verror ) );
%%PLOT
if plots
  figure()
  plot(x,postPredEval, '-r'); hold on;
  plot(x,tpdfEval, '-b')
  hold off
end 
