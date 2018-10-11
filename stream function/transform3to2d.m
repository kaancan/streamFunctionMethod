%% Transform 3d coil geometry to 2d
% 23.04.2018 - M. Kaan Can
% Maps the 3d coil geometry vertices and faces into 2d surface. Specific to
% helmet geometry, but may also work on similar shapes like a sphere.
% surface = 
%     x: X coordinates of coil, ncoilx1 array
%     y: Y coordinates of coil, ncoilx1  array
%     z: Z coordinates of coil, ncoilx1  array
%     n: normal vectors of each vertice, ncoilx3 matrix
%     faces: nfacesx3 matrix connecting vertices to form faces
%%
function surface2d = transform3to2d(surface)
    ncoil = length(surface.x);
    for i = 1:ncoil
        r(i) = sqrt(surface.x(i)^2 + surface.y(i)^2 +(surface.z(i)-max(surface.z))^2);
        ang(i) = atan2d(surface.y(i), surface.x(i));
        ang2(i) = atan2d(surface.z(i)-max(surface.z),sqrt(surface.x(i)^2 + surface.y(i)^2) );
    end
    surface2d.x = cosd(ang) .* r;
    surface2d.y = sind(ang) .* r;
    surface2d.faces = surface.faces;
    surface2d.centroids.x = mean(surface2d.x(surface2d.faces),2);
    surface2d.centroids.y = mean(surface2d.y(surface2d.faces),2);
    surface2d.ang = ang2;
end