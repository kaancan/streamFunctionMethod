function contourPlotsv2(surface, phi, n)
%%
surface2d = transform3to2d(surface);
[xq,yq] = meshgrid(-210:210,-235:122 );
vq = griddata(surface2d.x,surface2d.y,phi,xq,yq);
theta = griddata(surface2d.x,surface2d.y,surface2d.ang,xq,yq);
%%
mm(1) = min(phi); mm(2) = max(phi);
phi = phi - (mm(1)+mm(2))/2;
stp = (mm(2)-mm(1))/(n);

for i = 1:n
clear px py
l = mm(1) + stp/2 + (i-1) * stp;
c1 = contourc(vq, [l l]);

px = c1(1,:) - 211;
py = c1(2,:) - 235;

pr = sqrt(px.^2 +py.^2);
pr = pr';

if(sum(c1(1,:) == l) == 1)
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

    if(l < 0)
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
    parts = c1(2,(c1(1,:) == l));
    bin = zeros(length(parts), max(parts));
    ind = find(c1(1,:) == l);

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

        
        if(l < 0)
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
end