%% Discretize the current paths and visualize in both 2d and 3d
% 23.04.2018 - M. Kaan Can
% Calculates and plots the discrete current paths. Specific to helmet type 
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
%%
% function contourPlots(surface, In, n)
% In = edgeCorr * est_In; n = 3;
In = phi; n = 12;
surface2d = transform3to2d(surface);
[xq,yq] = meshgrid(-210:210,-235:122 );
vq = griddata(surface2d.x,surface2d.y,In,xq,yq);
theta = griddata(surface2d.x,surface2d.y,surface2d.ang,xq,yq);

%%
figure
[c] = contourc(vq,n); axis off; title('2D curent path distribution')
l = h.LevelList;
figure(50); hold on; axis off; view([150 30]); title('3D current path distribution')
for i = 1:length(l)
clear px py
c1 = contourc(vq, [l(i) l(i)]);

px = c1(1,:) - 211;
py = c1(2,:) - 235;
% 
% py(abs(px + 211) < 0.5) = [];
% px(abs(px + 211) < 0.5) = []; 
%  
pr = sqrt(px.^2 +py.^2);
pr = pr';


if(sum(c1(1,:) == l(i)) == 1)
    px(1) = [];
    py(1) = [];
    pr(1) = [];
    
    pang1 = atan2d(py,px);
    pang2 = diag(theta(round(py + 237),round(px + 211)));

    c_3d.x =  (px./abs(px))' .* sqrt((pr.*cosd(pang2)).^2 ./ (1+(tand(pang1)).^2)');
    c_3d.y = c_3d.x .* tand(pang1)';
    c_3d.z = pr .* sind(pang2) + max(surface.z);

    c_3d.x(pang1 == 90) = 0;
    c_3d.y(pang1 == 90) = sqrt(pr(pang1 == 90).^2 - (c_3d.z(pang1 == 90)-max(surface.z)).^2);
    c_3d.x(pang1 == -90) = 0;
    c_3d.y(pang1 == -90) = -sqrt(pr(pang1 == -90).^2 - (c_3d.z(pang1 == -90)-max(surface.z)).^2)';
% 
%     c_3d.x(1) = [];
%     c_3d.y(1) = [];
%     c_3d.z(1) = [];
    
    if(l(i) < 0)
        figure(50);
        line(c_3d.x,c_3d.y,c_3d.z,'Color',[0,0,1],'LineWidth',1)
        figure(2)
        line(px,py,'Color',[0,0,1],'LineWidth',1)
    else
        figure(50)
        line(c_3d.x,c_3d.y,c_3d.z,'Color',[1,0,0],'LineWidth',1)
        figure(2)
        line(px,py,'Color',[1,0,0],'LineWidth',1)
    end
    
else
    parts = c1(2,(c1(1,:) == l(i)));
    bin = zeros(length(parts), max(parts));
    ind = find(c1(1,:) == l(i));

    for j = 1:length(parts)
        bin(j,1:parts(j)) = ind(j)+1:ind(j)+parts(j);
        pang1 = atan2d(py(bin(j,find(bin(j,:)))),px(bin(j,find(bin(j,:)))));
        pang2 = diag(theta(round(py(bin(j,find(bin(j,:)))) + 236),round(px(bin(j,find(bin(j,:)))) + 211)));
        
        c_3d.x =  (px(bin(j,find(bin(j,:))))./abs(px(bin(j,find(bin(j,:))))))' .* sqrt((pr(bin(j,find(bin(j,:)))).*cosd(pang2)).^2 ./ (1+(tand(pang1)).^2)');
        c_3d.y = c_3d.x .* tand(pang1)';
        c_3d.z = pr(bin(j,find(bin(j,:)))) .* sind(pang2) + max(surface.z);
    
        c_3d.x(pang1 == 90) = 0;
        c_3d.y(pang1 == 90) = sqrt(pr(bin(j,pang1 == 90)).^2 - (c_3d.z(pang1 == 90)-max(surface.z)).^2);
        c_3d.x(pang1 == -90) = 0;
        c_3d.y(pang1 == -90) = -sqrt(pr(bin(j,pang1 == -90)).^2 - (c_3d.z(pang1 == -90)-max(surface.z)).^2)';

        
        if(l(i) < 0)
            figure(50);
            line(c_3d.x,c_3d.y,c_3d.z,'Color',[0,0,1],'LineWidth',1)
            figure(2)
            line(px(bin(j,find(bin(j,:)))),py(bin(j,find(bin(j,:)))),'Color',[0,0,1],'LineWidth',1)
        else
            figure(50);
            line(c_3d.x,c_3d.y,c_3d.z,'Color',[1,0,0],'LineWidth',1)
            figure(2)
            line(px(bin(j,find(bin(j,:)))),py(bin(j,find(bin(j,:)))),'Color',[1,0,0],'LineWidth',1)
        end
    end
end
end
% end