%% Interpolating helmet coil
% 01.03.2018 - M. Kaan Can
% Use 102-node helmet coil to interpolate a more detailed, 1467-node helmet
% shaped coil geometry.
% Updated on 30.07.2018
%%
close all;
clear all;
%% Load the coil and set up the variables
load helmet102.mat;
surface = struct('x',0,'y',0,'z',0,'n',[0;0;0;0;0]);
faces = FM_orig';

surface.x = VM_orig(:,1)*1000;surface.y = VM_orig(:,2)*1000; surface.z = VM_orig(:,3)*1000;
surface.n = normal_orig;

% Loop for number of iterations. 2 iterations give a coil with 1467 nodes, 2832 faces.
for j = 1:2
    tic
    nfaces = length(faces);
    ncoil = length(surface.x);
    % Find centre of each triangle
    cents.x = [mean([surface.x(faces(:,1)),surface.x(faces(:,2))],2);...
                mean([surface.x(faces(:,1)),surface.x(faces(:,3))],2);...
                mean([surface.x(faces(:,3)),surface.x(faces(:,2))],2)];
    cents.y = [mean([surface.y(faces(:,1)),surface.y(faces(:,2))],2);...
                mean([surface.y(faces(:,1)),surface.y(faces(:,3))],2);...
                mean([surface.y(faces(:,3)),surface.y(faces(:,2))],2)];
    cents.z = [mean([surface.z(faces(:,1)),surface.z(faces(:,2))],2);...
                mean([surface.z(faces(:,1)),surface.z(faces(:,3))],2);...
                mean([surface.z(faces(:,3)),surface.z(faces(:,2))],2)];
    
    % Find the normal vector at each centre by averaging 
    cents.n = [mean([surface.n(faces(:,1),1),surface.n(faces(:,2),1)],2), ...
        mean([surface.n(faces(:,1),2),surface.n(faces(:,2),2)],2), ...
        mean([surface.n(faces(:,1),3),surface.n(faces(:,2),3)],2) ;...
                mean([surface.n(faces(:,1),1),surface.n(faces(:,3),1)],2),...
        mean([surface.n(faces(:,1),2),surface.n(faces(:,3),2)],2),...
        mean([surface.n(faces(:,1),3),surface.n(faces(:,3),3)],2);...
                mean([surface.n(faces(:,3)),surface.n(faces(:,2))],2),...
        mean([surface.n(faces(:,3),2),surface.n(faces(:,2),2)],2),...
        mean([surface.n(faces(:,3),3),surface.n(faces(:,2),3)],2)];
    
    ncents = length(cents.x);
    
    % Create new coil nodes
    surfaceNew.x = [surface.x; cents.x];
    surfaceNew.y = [surface.y; cents.y];
    surfaceNew.z = [surface.z; cents.z];
    surfaceNew.n = [surface.n; cents.n];
    
    %Create new coil faces
    for i = 1:nfaces
        facesNew(4*(i-1)+1,:) = [faces(i,1),ncoil+i,ncoil+(1*ncents/3)+i];
        facesNew(4*(i-1)+2,:) = [ncoil+i,faces(i,2),ncoil+(2*ncents/3)+i];
        facesNew(4*(i-1)+3,:) = [ncoil+(2*ncents/3)+i,ncoil+(1*ncents/3)+i,faces(i,3)];
        facesNew(4*(i-1)+4,:) = [ncoil+i,ncoil+(2*ncents/3)+i,ncoil+(1*ncents/3)+i];
    end
    
    % Cleaning
    faces = facesNew; clear facesNew;
    surface = (surfaceNew); clear surfaceNew;
    nfaces = length(faces);
    ncoil = length(surface.x);
    
    % Eliminate the duplicate faces.  
    [n,  bin] = histc(surface.x, unique(surface.x));
    faces2 = zeros(size(faces(:,1)));
    for i = 1:length(faces)
        faces2(i,1) = find(bin == bin(find(surface.x == surface.x(faces(i,1)),1)),1);
        faces2(i,2) = find(bin == bin(find(surface.x == surface.x(faces(i,2)),1)),1);
        faces2(i,3) = find(bin == bin(find(surface.x == surface.x(faces(i,3)),1)),1);
    end
    faces3 = faces2(:);
    uni = unique(faces3);
    for k = 1:length(uni)
            faces3(faces3 == uni(k)) = k;
    end
    faces = reshape(faces3, nfaces,3);

    % Eliminate duplicate points
    surfNew = unique([surface.x, surface.y, surface.z],'rows','stable');
    surface.n = unique(surface.n,'rows','stable');
    surface.x = surfNew(:,1);
    surface.y = surfNew(:,2);
    surface.z = surfNew(:,3);
    surface.faces = faces;
    toc
end
%%
clearvars -except surface