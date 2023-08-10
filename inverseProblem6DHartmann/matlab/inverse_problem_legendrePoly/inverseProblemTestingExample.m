clear all; close all; clc;
addpath('../posteriorPredictiveApproach','-end')
%% POINTS
nres = 100;
xres = linspace(-1,1,nres);
%% MODEL
% Normal Density
sigma = 0.3.^2;
phi   = @(x) exp(-x.^2./2/sigma)/sqrt(2*pi*sigma); 
% Target Density
tpdf  = @(x) (phi(x-0.5)+phi(x+0.5))./2;          
% Normalized target density evaluations
tpdfOVERxres          = tpdf(xres);
normalizedTrueDensity = tpdfOVERxres./sum(tpdfOVERxres);
% Model approach functions evaluations
modelApproachFunctions = @(x) legendreApproachFunctions(x);
modelSigma2            = 0.01^2;
%% PARAMETERS
params = {
    {'theta1', 0.5, -10, 10}
    {'theta2', 0.5, -10, 10}
    {'theta3', 0.5, -10, 10}
    {'theta4', 0.5, -10, 10}
    {'theta5', 0.5, -10, 10}
    {'theta6', 0.5, -10, 10}
    {'theta7', 0.5, -10, 10}
    {'theta8', 0.5, -10, 10}
    {'theta9', 0.5, -10, 10}
    {'theta10', 0.5, -10, 10}
};
%% OPTIONS
options.method      = 'am';
options.nsimu       = 100000;
options.updatesigma = 0;
nord                = length(params);
% Proposal covariance
options.qcov        = eye(nord).*0.05^2 + (ones(nord)-eye(nord)).*0.01.^2;
burning             = 1000;
thin                = 10;
%% POSTERIOR PREDICTIVE APPROACH
predictedTheta = posteriorPredictiveApproach(xres, normalizedTrueDensity, modelApproachFunctions, modelSigma2, params, options, burning, thin, true);
