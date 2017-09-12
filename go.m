% this script performs slice-based dynamic shimming using one set of static
% second-order shims (xy, yz, xz, x2-y2, z2) and a separate set of dynamic shims
% (x, y, z, df) for each slice. it accounts for imperfections in the static
% second-order shims by a basis fitting and inversion process.

% compute coil calibration coefficients ........................................
 
% load data
%
% `o' is a cell array that contains the field maps (in order):
%   baseline (all shims = 0)
%   xy   = +3.5 A
%   zy   = +3.5 A
%   zx   = +3.5 A
%   x2y2 = +3.5 A
%   z2   = +3.5 A
% 
% `m' is a 4-D array that contains the corresponding magnitude images in the 4th
% dimension of the array
% 
% `c' contains the (x,y,z) center positions of each slice

clear, load coil

nshim = size(m,4); % # shims total
A = 3.5; % amps used for coil calibration

% get mask for each magnitude image
s = zeros(size(m));
for i = 1:nshim
  s(:,:,:,i) = mask2(m(:,:,:,i));
end
s = prod(s,4); % get intersection of masks for fitting spherical harmonics
s = imerode(s, strel('disk',4));

% mask field maps
for i = 1:nshim
  o{i} = s.*o{i};
end

% compute coil calibration matrix
M = coilm(o, s, A, c(:,3));

% perform slice-based dynamic shimming .........................................

clear, load M, load shim % load shim data

slshim(o1, s1, c1, 22, o2, s2, z2, 30, M, 'show', 1, 'out', 1);

