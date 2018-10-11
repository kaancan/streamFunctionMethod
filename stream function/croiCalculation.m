%% Stream function method part 3: Curl operator in ROI
% 28.02.2018 - M. Kaan Can
% The code that calculates the matrix for performing curl operation on the
% region of interest. The region is assumed to be a uniformly spaced
% rectangular prism. No extra functions are used.
% Inputs: 
% ROI geometry information:
% ROI = 
%     x: X coordinates of roi, nROI x 1 array
%     y: Y coordinates of roi, nROI x 1  array
%     z: Z coordinates of roi, nROI x 1  array
%     nsize: size of ROI in each dimension, 3 x 1 array
function Croi = croiCalculation(ROI)
tic
fprintf('Croi calculations...\n')
%% Croi calculations
nROIx = ROI.nsize(1);
nROIy = ROI.nsize(2);
nROIz = ROI.nsize(3);

px = eye(nROIx); 
nx = circshift(-px, [0, 1]);
Croinew.X = px+nx;
Croinew.X(end,:) = circshift(Croinew.X(end,:), [0, -1]);
Croinew.X = repmat({Croinew.X}, 1, nROIy*nROIz);
Croinew.X = blkdiag(Croinew.X{:});

py = eye(nROIx*nROIy); 
ny = circshift(-py, [0, nROIx]);
Croinew.Y = py+ny;
Croinew.Y(end-nROIx:end,:) = circshift(Croinew.Y(end-nROIx:end,:), [0, -nROIx]);
Croinew.Y = repmat({Croinew.Y}, 1, nROIz);
Croinew.Y = blkdiag(Croinew.Y{:});

Croi = [ -Croinew.Y, Croinew.X];
fprintf('Done.\n')
toc
end