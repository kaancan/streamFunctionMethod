%% Define a cylinder surface
% 23.04.2018 - M. Kaan Can
% A function to define a cylindrical surface, including the faces and
% normal vectors.
% Sizes are given in mm. Always check the coregistration with the ROI.
%%
function surface = defCoilSurface(cylinderRadius, cylinderLength, simStep, angStep)
surface = struct('x', 0, 'y', 0, 'z', 0, 'n', [0; 0; 0]);

theta = 0:angStep:359;
X = cylinderRadius*cosd(theta);
Y = cylinderRadius*sind(theta);
Z = 0:simStep:(cylinderLength)-simStep;

Z = repmat(Z, 1, length(X)); Z = sort(Z);   
X = repmat(X, 1, cylinderLength/(simStep));
Y = repmat(Y, 1, cylinderLength/(simStep));

surface.x = X(:); surface.y = Y(:); surface.z = Z(:);

N = zeros(3,1);
for i = 1:360/angStep
    N(:, i) = [X(i)/sqrt(X(i)^2 + Y(i)^2),Y(i)/sqrt(X(i)^2 + Y(i)^2),0];
end

N = repmat(N,1,cylinderLength/(simStep));
surface.n = N;
% surface.n(4,:) = atan2d(surface.n(2,:),surface.n(1,:));
surface.n = surface.n';
% quiver3(surface.x, surface.y, surface.z,surface.n(1,:)', surface.n(2,:)', surface.n(3,:)')  %Visualize n vectors

faces = delaunay(surface.x.*(surface.z+30), surface.y.*(surface.z+30));
bot = [];
for i = 1:length(faces)
    if(surface.z(faces(i,:)) == min(surface.z))
        bot = [bot i];
    end
end
faces(bot,:) = [];
surface.faces = faces;
end