%% Create a uniform cubic ROI
% 30.07.2018 - M. Kaan Can
% A function that defines the coordinates of a cubic ROI with given size.
% Use 'move' array to set the location.
%%
function ROI = defROI(sizeROI, stepROI, move)
if(nargin < 3)
    move = [0, 0, 0];
end
ROI = struct('x', 0, 'y', 0, 'z', 0);
ROI.x = repmat(-sizeROI/2:stepROI:sizeROI/2-stepROI, 1, (sizeROI/stepROI)^2)';
ROI.y = repmat(sort(repmat(-sizeROI/2:stepROI:sizeROI/2-stepROI, 1, sizeROI/stepROI)), 1, sizeROI/stepROI)';
ROI.z = sort(repmat(-sizeROI/2:stepROI:sizeROI/2-stepROI, 1, (sizeROI/stepROI)^2))';

ROI.x = ROI.x + move(1);
ROI.y = ROI.y + move(2);
ROI.z = ROI.z + move(3);

ROI.nsize = [sizeROI/stepROI sizeROI/stepROI sizeROI/stepROI];
end