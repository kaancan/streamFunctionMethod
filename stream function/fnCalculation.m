%% Stream function method part 1: fn Calculation
% 28.02.2018 - M. Kaan Can
% The code that calculates fn, which uses coil geometry to calculate the
% current density at each face of the coil when multiplied with the stream
% function. Uses the following other functions:
% detectEdge
% Inputs: 
% Coil geometry information:
% surface = 
%     x: X coordinates of coil, ncoil x 1 array
%     y: Y coordinates of coil, ncoil x 1 array
%     z: Z coordinates of coil, ncoil x 1 array
%     n: normal vectors of each vertice, ncoil x 3 matrix
%     faces: nfaces x 3 matrix connecting vertices to form faces
% Outputs:
% fn:           matrix that calculated current basis function in 2d
% fn3d:         matrix that calculated current basis function in 3d
% edgeCorr:     matrix that has the information related to edge nodes
function [fn,fn3d,edgeCorr] = fnCalculation(surface)
tic
fprintf('fn calculations...\n')
%% Calculate centroids and array lengths for future use
centroids.x = mean(surface.x(surface.faces),2);
centroids.y = mean(surface.y(surface.faces),2);
centroids.z = mean(surface.z(surface.faces),2);

ncoil = length(surface.x);
nfaces = length(surface.faces);
%% Edge Correction
% Detect edge nodes and create an edge correction matrix that enforces the
% stream function values to be equal on surface boundaries.
edgeNodes = detectEdgev2(surface); % See function for details.
noe = size(edgeNodes,1);
edgeCorr = zeros(ncoil,ncoil - length(edgeNodes(:)) + noe);
k = 0;
for i = 1:ncoil
    if(~sum(edgeNodes == i))
        edgeCorr(i,i-k) = 1;
    else
        for nn = 1:noe
            if(sum(edgeNodes(nn,:) == i))
                edgeCorr(i, ncoil - length(edgeNodes(:)) + nn) = 1;
                k = k + 1;
            end
        end
    end
end
%% Calculate fn 
neighbours = 0;
Ccoil = struct('all', 0, 'X', 0, 'Y', 0, 'Z', 0);
Ccoil.X = zeros(nfaces,ncoil);
Ccoil.Y = zeros(nfaces,ncoil);
Ccoil.Z = zeros(nfaces,ncoil);

for i = 1:nfaces
    ff = surface.faces(i,1);
    area = [];

    v1 = [surface.x(surface.faces(i,1))-surface.x(surface.faces(i,2)),surface.y(surface.faces(i,1))-surface.y(surface.faces(i,2)),surface.z(surface.faces(i,1))-surface.z(surface.faces(i,2))];
    v2 = [surface.x(surface.faces(i,1))-surface.x(surface.faces(i,3)),surface.y(surface.faces(i,1))-surface.y(surface.faces(i,3)),surface.z(surface.faces(i,1))-surface.z(surface.faces(i,3))];
    area = 0.5*norm(cross(v1,v2));
 
    n = [surface.n(surface.faces(i,1),1), surface.n(surface.faces(i,1),2), surface.n(surface.faces(i,1),3)];
    c1 = [centroids.x(i) - surface.x(surface.faces(i,1)), centroids.y(i) - surface.y(surface.faces(i,1)), centroids.z(i) - surface.z(surface.faces(i,1))];
    c2 = [surface.x(surface.faces(i,2)) - surface.x(surface.faces(i,1)), surface.y(surface.faces(i,2)) - surface.y(surface.faces(i,1)), surface.z(surface.faces(i,2)) - surface.z(surface.faces(i,1))];
    c3 = [surface.x(surface.faces(i,3)) - surface.x(surface.faces(i,1)), surface.y(surface.faces(i,3)) - surface.y(surface.faces(i,1)), surface.z(surface.faces(i,3)) - surface.z(surface.faces(i,1))];
    Ang(1) = dot(n,cross(c1,c2)); 
    Ang(2) = dot(n,cross(c1,c3));
    [a, b] = max(Ang);
    p = [2,3];
    p(b) = [];
    ff = [ff surface.faces(i, b+1) surface.faces(i, p)];
    
    pp = [surface.x(ff(3)) - surface.x(ff(2)), surface.y(ff(3)) - surface.y(ff(2)), surface.z(ff(3)) - surface.z(ff(2))];
    qq = [surface.x(ff(1)) - surface.x(ff(3)), surface.y(ff(1)) - surface.y(ff(3)), surface.z(ff(1)) - surface.z(ff(3))];
    ss = [surface.x(ff(2)) - surface.x(ff(1)), surface.y(ff(2)) - surface.y(ff(1)), surface.z(ff(2)) - surface.z(ff(1))];
    
    Ccoil.X(i, ff) = [pp(1), qq(1), ss(1)]/(2*area);
    Ccoil.Y(i, ff) = [pp(2), qq(2), ss(2)]/(2*area); 
    Ccoil.Z(i, ff) = [pp(3), qq(3), ss(3)]/(2*area);
end
%% Form the matrices
fn3d = [Ccoil.X; Ccoil.Y; Ccoil.Z] * edgeCorr;
fn = [Ccoil.X; Ccoil.Y] * edgeCorr;
fprintf('Done.\n')
toc
end