%% Stream function method part 2: Surface Integral over the coil faces
% 28.02.2018 (updated on 23.04.2018) - M. Kaan Can
% The code that calculates the matrix for taking the surface integral over
% the faces of coil. 
% Inputs:
% Surface geometry information:
% surface = 
%     x: X coordinates of coil, ncoil x 1 array
%     y: Y coordinates of coil, ncoil x 1  array
%     z: Z coordinates of coil, ncoil x 1  array
%     n: normal vectors of each vertice, ncoil x 3 matrix
%     faces: nfaces x 3 matrix connecting vertices to form faces
% ROI geometry information:
% ROI = 
%     x: X coordinates of roi, ncoil x 1 array
%     y: Y coordinates of roi, ncoil x 1  array
%     z: Z coordinates of roi, ncoil x 1  array
%     nsize: size of ROI in each dimension, 3 x 1 array
function [RR, RRxy] = surfintCalculations(surface,ROI)
tic
fprintf('Surface integral calculations...\n')
%% Calculate centroids and array lengths for future use
centroids.x = mean(surface.x(surface.faces),2);
centroids.y = mean(surface.y(surface.faces),2);
centroids.z = mean(surface.z(surface.faces),2);

ncoil = length(surface.x);
nfaces = length(surface.faces);
%% Calculating the matrix

const = 1;%1.2566*1e-6/(4*pi);
nROI = length(ROI.x);
for j = 1:nfaces
    v1 = [surface.x(surface.faces(j,1))-surface.x(surface.faces(j,2)),surface.y(surface.faces(j,1))-surface.y(surface.faces(j,2)),surface.z(surface.faces(j,1))-surface.z(surface.faces(j,2))];
    v2 = [surface.x(surface.faces(j,1))-surface.x(surface.faces(j,3)),surface.y(surface.faces(j,1))-surface.y(surface.faces(j,3)),surface.z(surface.faces(j,1))-surface.z(surface.faces(j,3))];
    area(j) = 0.5*norm(cross(v1,v2));
end
for i = 1:nROI
    for j = 1:nfaces
        R(i,j) = const * area(j)/(sqrt((centroids.x(j)-ROI.x(i))^2 + (centroids.y(j)-ROI.y(i))^2 + (centroids.z(j)-ROI.z(i))^2));
    end
end
%% Combining matrices
RR = blkdiag(R,R,R);
RRxy = blkdiag(R,R);
fprintf('Done.\n')
toc
end