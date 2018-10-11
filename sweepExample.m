clear all
close all
%% Load data
load('helmet1415.mat')              % Load Helmet geometry
load('headROI')                     % Load ROI geometry
load('symmetricTargetField1200')    % Load targetField
load('BrainMask1200.mat')
%% Use streamFunctionSweep function to perform a lambda sweep with desired constraint
[fitPercentage, residualInhomogeneity] = streamFunctionSweep(surface, ROI, BB, 'laplacian', [-5 5 2], 'log', w_brainMask)