%% Visualize the stream function and current distribution in both 2d and 3d
% 23.04.2018 - M. Kaan Can
% Plots the stream function distribution in 3d geometry. Maps the
% distribution to a 2D surface. Specific to helmet type coil.
% Inputs: 
% Coil geometry information:
% In: stream function, ncoilx1 array
% surface = 
%     x: X coordinates of coil, ncoilx1 array
%     y: Y coordinates of coil, ncoilx1  array
%     z: Z coordinates of coil, ncoilx1  array
%     n: normal vectors of each vertice, ncoilx3 matrix
%     faces: nfacesx3 matrix connecting vertices to form faces
% J: Current vectors. 3*nfacesx1 array. Optional.
% Updated on 31.07.2018: 
% 2d visualization is commented out for sake of generality.
%%
function visualization_In(In, surface, J)
    a = max(In);
    b = min(In);
    In = In - (a+b)/2;
    mm = max(abs(max(In)),abs(min(In)));
    figure
    hp = patch('faces',surface.faces,'vertices',[surface.x surface.y surface.z],'facevertexcdata',In/mm,'facecolor','interp','EdgeColor', 'interp');
    if(nargin == 3)
        nfaces = length(surface.faces);
        centroids.x = mean(surface.x(surface.faces),2);
        centroids.y = mean(surface.y(surface.faces),2);
        centroids.z = mean(surface.z(surface.faces),2);
        hold on
        quiver3(centroids.x, centroids.y, centroids.z, J(1:nfaces), J(nfaces+1: 2*nfaces), J(2*nfaces+1: 3*nfaces),'r')
    end
    axis off
    axis tight
    colormap(parula(12));colorbar;
    title('3D stream function and current distribution')
    
    
%     surface2d = transform3to2d(surface);
%     figure
%     hp = patch('faces',surface2d.faces,'vertices',[surface2d.x; surface2d.y]','facevertexcdata',In/mm,'facecolor','interp','EdgeColor', 'interp');
% %     caxis([-1 1])
%     axis off
%     axis tight
%     colormap(parula);%colorbar;
%     title('2D stream function distribution')
end

