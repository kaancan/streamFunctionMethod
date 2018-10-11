clear all
close all
%% Load data
load('helmet1415.mat')              % Load coil geometry
load('headROI')                     % Load ROI geometry
load('symmetricTargetField1200')    % Load targetField
load('BrainMask1200.mat')           % Load mask for ROI
%% Use streamFunction function to obtain a phi distribution
[phi, estField] = streamFunction(surface, ROI, BB, 'laplacian', 1e3,w_brainMask);
%% Evaluation
fieldImprovementRatio = 100*(1-sum(w_brainMask.*(estField - BB).^2)/ sum((BB).^2))
standardDeviation = 100*(1-std(w_brainMask.*(estField-BB))/std(w_brainMask.*BB))
[residuals, originals] = resInhomogeneityCalculation(estField, BB)