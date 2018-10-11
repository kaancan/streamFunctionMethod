%% Discretize the current paths and visualize in both 2d and 3d
% 23.04.2018 - M. Kaan Can
% Calculates and plots the discrete current paths. Specific to cylindrical 
% coil.
% Inputs: 
% In: stream function, ncoilx1 array
% surface = 
%     x: X coordinates of coil, ncoilx1 array
%     y: Y coordinates of coil, ncoilx1  array
%     z: Z coordinates of coil, ncoilx1  array
%     n: normal vectors of each vertice, ncoilx3 matrix
%     faces: nfacesx3 matrix connecting vertices to form faces
% n: number of contours, scalar
%% Create 2d contour map
ncoil = length(surface.x);
surface2d = transform3to2cylinder(surface);
nedgenodes = 144;%sum(sum(edgeCorr(:,end) == 1));
n = 12;
In = reshape(phi,nedgenodes,1872/nedgenodes); 
In = In';
% % Duplicate the boundary for continuity
In = [In In(:,1)];
[xq,yq] = meshgrid(min(surface2d.ang):max(surface2d.ang),min(surface2d.h):max(surface2d.h));
vq = griddata(surface2d.ang,surface2d.h,phi,xq,yq);
vq = [vq vq(:,1)];
%%
mm(1) = min(phi); mm(2) = max(phi);
phi = phi - (mm(1)+mm(2))/2;
stp = (mm(2)-mm(1))/(n);

figure; 
[c,h] = contour(vq,12);
l = h.LevelList;
figure(50); hold on;
coil = struct('x',[],'y',[],'z',[]);
J = struct('x',[],'y',[],'z',[]);
for i = 1:n
    
l(i) = mm(1) + stp/2 + (i-1) * stp;
c1 = contourc(vq, [l(i) l(i)]);
% r = 8;
c_ang = (146/nedgenodes) .* c1(1,:);

c_3d.x = cosd(c_ang);
c_3d.y = sind(c_ang);
c_3d.z = c1(2,:);

parts = c1(2,c1(1,:) == l(i));
bin = zeros(length(parts), max(parts));
ind = find(c1(1,:) == l(i));

    for j = 1:length(parts)
        bin(j,1:parts(j)) = ind(j)+1:ind(j)+parts(j);
%         figure(50)
        coil.x = [coil.x c_3d.x(bin(j,find(bin(j,:))))]; coil.x(end) = [];
        coil.y = [coil.y c_3d.y(bin(j,find(bin(j,:))))]; coil.y(end) = [];
        coil.z = [coil.z c_3d.z(bin(j,find(bin(j,:))))]; coil.z(end) = [];
        

        
        if(l(i) < l(6))
            figure(50);line(150*c_3d.x(bin(j,find(bin(j,:)))),150*c_3d.y(bin(j,find(bin(j,:)))),10*(-12.5 + c_3d.z(bin(j,find(bin(j,:))))),'Color',[0,0,1],'LineWidth',1)   
            figure(2); line(c1(1,bin(j,find(bin(j,:)))),c1(2,bin(j,find(bin(j,:)))),'Color',[0,0,1],'LineWidth',1)
%             J.x = [J.x -gradient(c_3d.x(bin(j,find(bin(j,:)))))]; J.x(end) = [];
%             J.y = [J.y -gradient(c_3d.y(bin(j,find(bin(j,:)))))]; J.y(end) = [];
%             J.z = [J.z -gradient(c_3d.z(bin(j,find(bin(j,:)))))]; J.z(end) = [];
        else
%             J.x = [J.x gradient(c_3d.x(bin(j,find(bin(j,:)))))]; J.x(end) = [];
%             J.y = [J.y gradient(c_3d.y(bin(j,find(bin(j,:)))))]; J.y(end) = [];
%             J.z = [J.z gradient(c_3d.z(bin(j,find(bin(j,:)))))]; J.z(end) = [];
            figure(50);line(150*c_3d.x(bin(j,find(bin(j,:)))),150*c_3d.y(bin(j,find(bin(j,:)))),10*(-12.5 + c_3d.z(bin(j,find(bin(j,:))))),'Color',[1,0,0],'LineWidth',1)             
            figure(2); line(c1(1,bin(j,find(bin(j,:)))),c1(2,bin(j,find(bin(j,:)))),'Color',[1,0,0],'LineWidth',1)
% 
        end
    end

end
%     scatter3(coil.x, coil.y, coil.z,'.')
%     figure
%     quiver3(coil.x,coil.y,coil.z,J.x,J.y,J.z)