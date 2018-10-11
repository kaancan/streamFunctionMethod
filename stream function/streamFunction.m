%% Stream function method
% 31.07.2018 - M. Kaan Can
% The code that calculates and visualizes the stream function distribution
% for a given surface, ROI and target magnetic field with a selected
% constraint.
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
% lambda = 1 x 1 numeric. Weighting factor for the selected constraint.
% W
% w = nROI x 1 vector. Weighting vector for the ROI. If not given, it will
% be taken 1.
% Outputs:
% Phi
% phi = ncoil x 1 vector. The resulting stream function distribution.
% Estimated Field
% estField = nROI x 1 vector. An estimation of magnetic field that is
% created by the resulting stream function. Substract from the target field
% in order to find the resulting field.
function [phi,estField] = streamFunction(surface, ROI, targetField, constraint, lambda, w)


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

if(~isnumeric(lambda) || ~isequal(size(lambda), size(1)))
    error('Invalid lambda value.')
end

if(nargin == 6)
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
[const, grad] = laplacian(surface,edgeCorr);
% Add new constraint function here!
% const = ;
%% Solve
wcn = cn.*repmat(w, [1 size(cn,2)]);
mat = (wcn' * wcn) + lambda*(const' * const);
mat_pinv = pinv(mat); 
rhs = (wcn' * targetField);
%% Set outputs
phi = edgeCorr * mat_pinv * rhs;
estField = wcn * mat_pinv * rhs;
%% Visualize results
visualization_In(phi, surface, fn3d * mat_pinv * rhs);
end