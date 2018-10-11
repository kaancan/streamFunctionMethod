%% Transform 3d coil geometry to 2d
% 23.04.2018 - M. Kaan Can
% Maps the 3d coil geometry vertices and faces into 2d surface. Specific to
% cylinderical coil.
% surface = 
%     x: X coordinates of coil, ncoilx1 array
%     y: Y coordinates of coil, ncoilx1  array
%     z: Z coordinates of coil, ncoilx1  array
%     n: normal vectors of each vertice, ncoilx3 matrix
%     faces: nfacesx3 matrix connecting vertices to form faces
%%
function surface2d = transform3to2cylinder(surface)
    surface2d.ang = atan2d(surface.y, surface.x)';
    surface2d.h = surface.z';
    surface2d.faces = surface.faces;
    nfaces = length(surface.faces);
%     for j = nfaces
%         if(sum(surface2d.ang(surface.faces(j,:)) == min(surface2d.ang)))% && sum(surface.faces(j,:) == max(surface2d.ang)))
% %             surface2d.faces(j,:) = [];
%             j
%         end
%     end
    surface2d.ang(surface2d.ang < -90) = surface2d.ang(surface2d.ang < -90) + 360;
    surface2d.faces(logical(sum(surface2d.ang(surface.faces) == min(surface2d.ang),2)) & logical(sum(surface2d.ang(surface.faces) == max(surface2d.ang),2)),:) = [];

end