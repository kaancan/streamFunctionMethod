%% Symmetric coil calculation
% 23.04.2018 - M. Kaan Can
% Interpolate and create a symmetric helmet coil using 102-node helmet
% coil. Needs "head_interpolation_KaanCan.m" file to work. The result is a
% symemtric and 1415-node helmet shaped coil geometry including faces and
% normals of the surface. Specific to this coil geometry.
%% 
clear all
close all
%% Load coil and interpret the detailed coil geometry from 102 node helmet
head_interpolation_KaanCan;
nfaces = length(surface.faces);
ncoil  = length(surface.x);
noe = 1;
surface.x = surface.x + 3;
surface.y = surface.y - 30;
surface.z = surface.z - 40;

surfaceSym.x = surface.x(surface.x > 0);
surfaceSym.y = surface.y(surface.x > 0);
surfaceSym.z = surface.z(surface.x > 0);
surfaceSym.n = surface.n(surface.x > 0,:);

surfaceSym.x((surfaceSym.x) < 9 & surfaceSym.x > 0) = 0;
surfaceSym.n((surfaceSym.x) < 9 & surfaceSym.x > 0,1) = 0;

k = 1;
for i = 1:nfaces
   if(sum(surface.x(surface.faces(i,:)) > 0) == 3) 
        surfaceSym.faces(k,:) = surface.faces(i,:);
        k = k+1;
   end
end
l = 0;
for i = 1:ncoil
    if(surface.x(i) < 0)
        l = l+1;
    end
    surfaceSym.faces(surfaceSym.faces == i) = surfaceSym.faces(surfaceSym.faces == i) - l;
end
%%
nn = length(surfaceSym.x);
nc = length(surfaceSym.faces);
surfaceSym.x(nn+1:2*nn) = -surfaceSym.x(1:nn);
surfaceSym.y(nn+1:2*nn) = surfaceSym.y(1:nn);
surfaceSym.z(nn+1:2*nn) = surfaceSym.z(1:nn);
surfaceSym.n(nn+1:2*nn,1) = -surfaceSym.n(1:nn,1); 
surfaceSym.n(nn+1:2*nn,2:3) = surfaceSym.n(1:nn,2:3); 

surfaceSym.faces(nc+1:2*nc,:) = surfaceSym.faces(1:nc,:) + nn;
%%
ind = find(surfaceSym.x == 0);
ind = ind(ind > nn);

for i = 1:length(ind)
    surfaceSym.faces(surfaceSym.faces == ind(i)) = surfaceSym.faces(surfaceSym.faces == ind(i)) - nn;

end

%%
surfaceSym.x(ind) = [];
surfaceSym.y(ind) = [];
surfaceSym.z(ind) = [];
surfaceSym.n(ind,:) = [];
k = 0;
for i = 1:1464
    if(sum(ind == i))
        surfaceSym.faces(surfaceSym.faces >= i-k) = surfaceSym.faces(surfaceSym.faces >= i-k)-1;
        k = k+1;
    end
end

%%
clear surface
clear faces
surface = surfaceSym;
%%
for i = 1:length(surface.x)
    norms = sqrt(surface.n(i,1).^2 + surface.n(i,2).^2 + surface.n(i,3).^2);
    surface.n(i,:) = surface.n(i,:)/norms;
end
%%
ncoil = length(surface.x);
nfaces = length(surface.faces);
clearvars -except surface