%% Stream function lambda sweep
% 31.07.2018 - M. Kaan Can
% The code that sweeps the weighting of the constraint in given logarithmic
% or linear range. Magnetic field fit percentage and residual
% inhomogeneities are returned as vectors.
% Inputs: 
% Coil geometry information:
% surface = 
%     x: X coordinates of coil, ncoil x 1 array
%     y: Y coordinates of coil, ncoil x 1 array
%     z: Z coordinates of coil, ncoil x 1 array
%     n: normal vectors of each vertice, ncoil x 3 matrix
%     faces: nfaces x 3 matrix connecting vertices to form faces
% ROI geometry information:
% ROI = 
%     x: X coordinates of roi, nROI x 1 array
%     y: Y coordinates of roi, nROI x 1  array
%     z: Z coordinates of roi, nROI x 1  array
%     nsize: size of ROI in each dimension, 3 x 1 array
% Target Field
% targetField = nROI x 1 array. If the target field is in 3D notation, can 
% be converted to proper notation by (:) function.
% Constraint
% constraint = a string field. !!! Only Laplacian constraint is implemented
% for now.!!! Any new constraint should be prepared as a script that
% calculates a related matrix operator. 
% Lambda
% lambda = 3 x 1 vector. Weighting factor for the selected constraint.
% first two elements of the vector determines the minimum and maximum
% limits of the sweep. If the sweep type is selected as logarithmic, the
% numbers corresponds to powers of 10. Last element is the number of points
% to be calculated in the range.
% Sweep Type
% sweepType = String. 'log' or 'lin' can be selected for logarithmic or
% linear sweep in the chosen range. 
% W
% w = nROI x 1 vector. Weighting vector for the ROI. If not given, it will
% be considered 1.
% Outputs
% Magnetic Fit Percentage
% fitPerc = nrange x 1 vector. Percentage value of fitting of target and
% resulting magnetic fields. Each value is the fit of the corresponding
% lambda value.
% Residual Inhomogeneity
% residuals = 
%       perc95:(nrange + 1) x 1 vector. Maximum residual inhomogeneity that 95% of the ROI has.
%       perc90:(nrange + 1) x 1 vector. Maximum residual inhomogeneity that 90% of the ROI has.
%       perc80:(nrange + 1) x 1 vector. Maximum residual inhomogeneity that 80% of the ROI has.
% Last value is corresponds to unshimmed target field. 
function [fitPerc, residuals] = streamFunctionSweep(surface, ROI, targetField, constraint, lambda, sweepType, w)

% Check the inputs
if(isstruct(surface))
    if(~isfield(surface,'x'))
        error('Surface does not have x component.')
    elseif(~isfield(surface,'y'))
        error('Surface does not have y component.')
    elseif(~isfield(surface,'z'))
        error('Surface does not have z component.')
    elseif(length(surface.x) ~= length(surface.y) || length(surface.x) ~= length(surface.z) || length(surface.z) ~= length(surface.y))
        error('Unmatched surface dimensions.')
    elseif(~isfield(surface,'n'))
        error('Surface does not have normal component.')
    elseif(~isfield(surface,'faces'))
        error('Surface does not have faces component.')
    end
else
    error('Surface input is not a struct.')
end

if(isstruct(ROI))
    if(~isfield(ROI,'x'))
        error('ROI does not have x component.')
    elseif(~isfield(ROI,'y'))
        error('ROI does not have y component.')
    elseif(~isfield(ROI,'z'))
        error('ROI does not have z component.')
    elseif(length(ROI.x) ~= length(ROI.y) || length(ROI.x) ~= length(ROI.z) || length(ROI.z) ~= length(ROI.y))
        error('Unmatched ROI dimensions.')
    elseif(~isfield(ROI,'nsize') || ~isequal(size(ROI.nsize), size([1 1 1])))
        error('ROI does not have nsize component or invalid nsize.')
    end
else
    error('ROI input is not a struct.')
end

if(isnumeric(targetField))
    if(length(targetField) ~= length(ROI.x))
        error('Invalid target field size.')
    end
else
    error('Target field should be numeric.')
end

if(~exist(constraint))
    error('.m file for the constraint is not found.')
end

if(~isnumeric(lambda) || ~isequal(size(lambda), size([1 1 1])))
    error('Invalid lambda value.')
end

if(sweepType == 'log')
    range = logspace(lambda(1), lambda(2), lambda(3));
elseif(sweepType == 'lin')
    range = linspace(lambda(1), lambda(2), lambda(3));
else
    error('Invalid sweep type.')
end

if(nargin == 7)
    if(~isnumeric(w) || ~isequal(length(ROI.x), length(w)))
        error('Invalid w.')
    end
else
    w = ones(size(ROI.x));
end
%% Definitions
ncoil   = length(surface.x);
nfaces  = length(surface.faces);
nROI    = length(ROI.x);
phi = zeros(ncoil,1);

%% Calculate magnetic field
[fn,fn3d,edgeCorr] = fnCalculation(surface);   % Calculate fn
[RR, RRxy] = surfintCalculations(surface,ROI); % Calculate surface integral on coil
Croi = croiCalculation(ROI);                   % Calculate curl operator in ROI
cn = (Croi * RRxy * fn);                       % Combine all matrices
%% Calculate constraint
[laplap, grad] = laplacian(surface,edgeCorr);
%%
wcn = cn.*repmat(w, [1 size(cn,2)]);
for i = 1:lambda(3)
    lamb = range(i);
    mat = (wcn' * wcn) + lamb*(laplap' * laplap);
    mat_pinv = pinv(mat); 
    rhs = (wcn' * targetField);
    estField = wcn * mat_pinv * rhs;
    fitPerc(i) = 100*(1-sum((estField - targetField).^2)/ sum((targetField).^2));
    [res, originals] = resInhomogeneityCalculation(estField, targetField);
    residuals.perc95(i) = res(1); residuals.perc95(lambda(3)+1) = originals(1);
    residuals.perc90(i) = res(2); residuals.perc90(lambda(3)+1) = originals(2);
    residuals.perc80(i) = res(3); residuals.perc80(lambda(3)+1) = originals(3);
    pp = 100*i/lambda(3);
    fprintf('%.2f percent completed.\n', pp)
end
end