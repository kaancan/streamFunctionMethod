%% Symmetric coil calculation
% 23.04.2018 - M. Kaan Can
% Loads and interpolates the off-resonance magnetic field. Resolution can
% be kept same or reduced(Comment out the other one). Also co-registers 
% with the coil.
%%
% clear all;
% close all;
%%
mri=load_untouch_nii('Mode1.nii');
vol=mri.img;
vol=permute(vol,[3 1 2]);
vol=flipdim(vol,3);
[x,y,z]=meshgrid([-38:37].*2, [-46:46].*2,  [-38:37].*2);
x=x+1;
y=y-20;
z=z+8;

%% Full resolution symmetry
% x0 = x(:,1:38,:);   x1 = x(:,39:76,:);
% y0 = y(:,1:38,:);   y1 = y(:,39:76,:);
% z0 = z(:,1:38,:);   z1 = z(:,39:76,:);
% v0 = vol(:,1:38,:); v1 = vol(:,39:76,:);
% 
% x0 = [x0, -flipdim(x0,2)];      x1 = [-flipdim(x1,2), x1];
% y0 = [y0,  y0];                 y1 = [ y1, y1];
% z0 = [z0,  z0];                 z1 = [ z1, z1];
% v0 = [v0,  flipdim(v0,2)];      v1 = [ flipdim(v1,2), v1];
% 
% vv = (v0 + v1)/2;
% vv = permute(vv,[2,1,3]);
% BB = vv(:);
% 
% ROI.size = [size(unique(x),1) size(unique(y),1) size(unique(z),1)];
% ROI.x = repmat(unique(x(:)),ROI.size(2)*ROI.size(3),1);
% ROI.y = repmat(sort(repmat(unique(y(:)), ROI.size(1), 1)), ROI.size(3), 1);
% ROI.z = z(:);
% nROI = length(ROI.x);
%% Lower the resolution symmetry
n = 8;
xx = x(n/4:n:end,n/4:n:end,n/4:n:end)+1;
yy = y(n/4:n:end,n/4:n:end,n/4:n:end);
zz = z(n/4:n:end,n/4:n:end,n/4:n:end);
vv = vol(n/4:n:end,n/4:n:end,n/4:n:end);

x0 = xx(:,1:5,:);   x1 = xx(:,6:10,:);
y0 = yy(:,1:5,:);   y1 = yy(:,6:10,:);
z0 = zz(:,1:5,:);   z1 = zz(:,6:10,:);
v0 = vv(:,1:5,:);   v1 = vv(:,6:10,:);

x0 = [x0, -flipdim(x0,2)];      x1 = [-flipdim(x1,2), x1];
y0 = [y0,  y0];                 y1 = [ y1, y1];
z0 = [z0,  z0];                 z1 = [ z1, z1];
v0 = [v0,  flipdim(v0,2)];      v1 = [ flipdim(v1,2), v1];

vv = (v0 + v1)/2;
%% Convert the notation to one can be used by stream function code.
vv = permute(vv,[2,1,3]);
BB = vv(:);

ROI.nsize = [size(unique(xx),1) size(unique(yy),1) size(unique(zz),1)];
ROI.x = repmat(unique(xx(:)),ROI.nsize(2)*ROI.nsize(3),1);
ROI.y = repmat(sort(repmat(unique(yy(:)), ROI.nsize(1), 1)), ROI.nsize(3), 1);
ROI.z = zz(:);
nROI = length(ROI.x);
%% 
clearvars -except ROI BB