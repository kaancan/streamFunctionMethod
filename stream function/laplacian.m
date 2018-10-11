%% Laplacian calculation
% 23.04.2018 - M. Kaan Can
% Function to calculate the laplacian operator. Also calculates gradient
% operator.
% Inputs: 
% Coil geometry information:
% surface = 
%     x: X coordinates of coil, ncoil x 1 array
%     y: Y coordinates of coil, ncoil x 1  array
%     z: Z coordinates of coil, ncoil x 1  array
%     n: normal vectors of each vertice, ncoil x 3 matrix
%     faces: nfaces x 3 matrix connecting vertices to form faces
% Outputs:
% laplap:   matrix that calculates the laplacian of the stream function over
% the coil surface.
% grad:     matrix that calculates the gradient of the stream function over 
% the coil surface.  
function [laplap, grad] = laplacian(surface,edgeCorr)
%% Calculate centroids and array lengths for future use
centroids.x = mean(surface.x(surface.faces),2);
centroids.y = mean(surface.y(surface.faces),2);
centroids.z = mean(surface.z(surface.faces),2);

ncoil = length(surface.x);
nfaces = length(surface.faces);
%% Partial derivatives
fprintf('Partial derivatives...\n')
tic
for i = 1:nfaces
    facenormals.x(i) = mean(surface.n(surface.faces(i,:),1),1);
    facenormals.y(i) = mean(surface.n(surface.faces(i,:),2),1);
    facenormals.z(i) = mean(surface.n(surface.faces(i,:),3),1);
end

facenormals.x = facenormals.x./sqrt((facenormals.x).^2 + (facenormals.y).^2 + (facenormals.z).^2);
facenormals.y = facenormals.y./sqrt((facenormals.x).^2 + (facenormals.y).^2 + (facenormals.z).^2);
facenormals.z = facenormals.z./sqrt((facenormals.x).^2 + (facenormals.y).^2 + (facenormals.z).^2);

facenormals.x = facenormals.x';
facenormals.y = facenormals.y';
facenormals.z = facenormals.z';

tr = triangulation(surface.faces, surface.x, surface.y, surface.z);
delta = 1e-10;

Pz0 = centroids.z - delta;
Pz1 = centroids.z + delta;
Py0 = centroids.y - delta;
Py1 = centroids.y + delta;
Px0 = centroids.x - delta;
Px1 = centroids.x + delta;

dz = (cartesianToBarycentric(tr, (1:nfaces)',[centroids.x centroids.y Pz1]) - cartesianToBarycentric(tr, (1:nfaces)',[centroids.x centroids.y Pz0]))/(2*delta);
dy = (cartesianToBarycentric(tr, (1:nfaces)',[centroids.x Py1 centroids.z]) - cartesianToBarycentric(tr, (1:nfaces)',[centroids.x Py0 centroids.z]))/(2*delta);
dx = (cartesianToBarycentric(tr, (1:nfaces)',[Px1 centroids.y centroids.z]) - cartesianToBarycentric(tr, (1:nfaces)',[Px0 centroids.y centroids.z]))/(2*delta);

dx(dx == 0) = delta;
dy(dy == 0) = delta;
dz(dz == 0) = delta;

Ccoil = struct('all', 0, 'X', 0, 'Y', 0, 'Z', 0);
for i = 1:nfaces
    Ccoil.X(i,surface.faces(i,:)) = dx(i,:);
    Ccoil.Y(i,surface.faces(i,:)) = dy(i,:);
    Ccoil.Z(i,surface.faces(i,:)) = dz(i,:);
end

Ccoil.alls = [Ccoil.X; Ccoil.Y; Ccoil.Z];
Ccoil.alll = [Ccoil.alls, zeros(3*nfaces,ncoil), zeros(3*nfaces,ncoil); zeros(3*nfaces,ncoil), Ccoil.alls, zeros(3*nfaces,ncoil); zeros(3*nfaces,ncoil), zeros(3*nfaces,ncoil), Ccoil.alls];
fprintf('Done.\n')
toc

%% Normal Correction - faces
fprintf('Normal correction...\n')
tic
nnfx = zeros(nfaces, 3*nfaces);

for i = 1:nfaces
    nnfx(i,i)           = facenormals.x(i);
    nnfx(i,i+nfaces)    = facenormals.y(i);
    nnfx(i,i+2*nfaces)  = facenormals.z(i);
end

nnf = zeros(3*nfaces, nfaces);

nnf(1:nfaces,:)           = diag(facenormals.x);
nnf(nfaces+1:2*nfaces,:)   = diag(facenormals.y);
nnf(2*nfaces+1:3*nfaces,:) = diag(facenormals.z);
subNormss = eye(3*nfaces) - nnf * nnfx;
subNormsl = blkdiag(subNormss,subNormss,subNormss);
fprintf('Done.\n')
toc

%% Face to vert
fprintf('face2ver calculation...\n')
tic
face2vert = zeros(3*ncoil,3*nfaces);
for i = 1:ncoil
    cnct = 0;
    
    for j = 1:nfaces
       if(sum(surface.faces(j,:) == i))
          cnct = [cnct j];
       end
    end
    denom = length(cnct) - 1;
    face2vert(i,cnct(2:end)) = 1/denom;
    face2vert(ncoil + i,nfaces+cnct(2:end)) = 1/denom;
    face2vert(2*ncoil + i,2*nfaces+cnct(2:end)) = 1/denom;
end
face2vert1d = face2vert(1:ncoil,1:nfaces);
fprintf('Done.\n')
toc
%% Select and sum matrices
fprintf('Select and summation matrices...\n')
tic
zzz = zeros(nfaces);
ooo = eye(nfaces);
slctLap = [ooo,zzz,zzz,zzz,zzz,zzz,zzz,zzz,zzz; ...
           zzz,zzz,zzz,zzz,ooo,zzz,zzz,zzz,zzz; ...
           zzz,zzz,zzz,zzz,zzz,zzz,zzz,zzz,ooo ];
sumsum = [ooo,ooo,ooo];       
fprintf('Done.\n')
toc
%% Combine matrices
fprintf('Combine matrices...\n')
tic
grad = face2vert * subNormss * Ccoil.alls  * edgeCorr;
laplap = sumsum * slctLap * subNormsl * Ccoil.alll * grad;
fprintf('Done.\n')
toc
end