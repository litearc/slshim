function [sShims dShims] = slshim(axB0, axMask, axXYZ, axFov, sgB0, sgMask, ...
  sgXYZ, sgFov, coilM, varargin)
  %
  %  minimizes the off-resonance in each slice in the least-squares sense.
  %  uses the  same second-order gradients for all slices, and varies the
  %  the linear gradients and frequency offset for each slice.
  %
  %  gradients:
  %    quadratic : xy, xz, yz, x2-y2, z2-.5(x2+y2)
  %    linear    : x, y, z
  %    constant  : bo
  %
  %  function [sShims dShims] = slshim(axB0, axMask, axXYZ, axFov, sgB0, ...
  %    sgMask, sgXYZ, sgFov, coilM, varargin)
  %
  %  inputs ....................................................................
  %  axB0             axial B0 map. [x y z] (Hz)
  %  axMask           mask for `oa'. [x y z] (binary)
  %  axXYZ            center position (x,y,z) for each slice in `oa'.
  %                   [slices (x,y,z)] (cm)
  %                   the z-locations are centered around 0 due to the table
  %                   shifting
  %  axFov            field of view for `oa'. (cm)
  %  sgB0             sagittal B0 map, 1 slice. [z x] (Hz)
  %  sgMask           mask for `os'. [z x]
  %  sgXYZ            position (x,y,z) of slice in `os' (cm)
  %  sgFov            field of view for 'os' (cm)
  %  coilM            matrix where M(i,j) is amount of field i that coil j
  %                   generates, e.g. M(1,2) is the amount of xy field generated
  %                   by the zy shim. M is [9 x 5], i.e. the first 5 columns
  %                   below, corresponding to the static shims, but is
  %                   automatically augmented to include the dynamic shims.
  % 
  %                   units of entries of M (augmented):
  % 
  %                        xy zy zx XY z2    x y z     bo
  %                      .----------------;--------;--------.
  %                   xy |                |        |        |
  %                   zy |     Hz/cm2     | Hz/cm2 | Hz/cm2 |
  %                   zx |     ------     | ------ | ------ |
  %                   XY |        A       |  G/cm  |   Hz   |
  %                   z2 |                |        |        |
  %                      |----------------+--------+--------|
  %                    x |      Hz/cm     |  Hz/cm | Hz/cm  |
  %                    y |      -----     |  ----- | -----  |
  %                    z |        A       |   G/cm |   Hz   |
  %                      |----------------+--------+--------|
  %                      |        Hz      |   Hz   |   Hz   |
  %                   bo |       ---      |  ----  |   --   |
  %                      |        A       |  G/cm  |   Hz   |
  %                      `----------------'--------'--------'
  %
  %  options ...................................................................
  %  writeFile        output shim file? (default = 0)
  %  gam              gamma (Hz/G) (default = 4258)
  %  objNorm          norm to use for objective function
  %                   options: 'L1', 'L2' (default)
  %  magWt            mag image used as weights in optimization
  %  pfitRangeZ       thickness in z (cm) used to do polynomial fit for
  %                   z-gradient calculation (default = 8)
  %  polyOrd          polynomial order for z-gradient calculation (default = 2)
  %  show             display various images (default = 0)
  %
  %  outputs ...................................................................
  %  sShims           second-order amplitudes: [xy zy zx x2-y2 z2-.5(x2+y2)]
  %  dShims           linear + offset for each slice: [x y z o] (nslices x 3)
  %
 
  % set default arguments
  v = ap2s(varargin);
  writeFile  = def(v, 'writeFile', 0);
  gam        = def(v, 'gam', 4258);
  objNorm    = def(v, 'objNorm', 'L2');
  magWt      = def(v, 'magWt', []);
  pfitRangeZ = def(v, 'pfitRangeZ', 8);
  polyOrd    = def(v, 'polyOrd', 2);
  show       = def(v, 'show', 0);

  % ............................................................................
  [ny, nx, nz] = size(axB0);
  nsShims = 3; % # static shims (excluding xz and yz)
  ndShims = 3; % # dynamic shims (excluding z, which is handled later)
  axB0 = axB0.*axMask; % apply masks
  sgB0 = sgB0.*sgMask;

  [xGrid,yGrid,zGrid] = get_xyz(axFov, axFov, nx, ny, axXYZ);
 
  % for optimization, we only care about the field inside mask
  idxMask = find(axMask);
  lMask = length(idxMask);

  % (x,y,z) contain the x, y, and z positions for each voxel within axial mask
  x = xGrid(idxMask);
  y = yGrid(idxMask);
  z = zGrid(idxMask);

  % 'in' contains the index of the first pixel inside the mask for each slice
  idx1stVoxInSlc = zeros(1,nz);
  for i = 1:nz
    idx1stVoxInSlc(i) = find(z==axXYZ(i,3),1);
  end

  % position encoding matrices for static and dynamic shims ....................
  Estatic = [x.*y, x.^2-y.^2, z.^2-(x.^2+y.^2)/2]; % don't include xz and yz
  E = [Estatic, zeros(lMask,ndShims*nz)];
  % for each slice, fill in one "block" section in the matrix
  % c1,c2 are the beginning/end columns
  % r1,r2 are the beginning/end rows
  for i = 1:nz
    c1 = nsShims+1+ndShims*(i-1);
    c2 = nsShims+ndShims*i;
    r1 = idx1stVoxInSlc(i);
    if i == nz, r2 = lMask; else r2 = idx1stVoxInSlc(i+1)-1; end
    E(r1:r2, c1:c2) = [x(r1:r2) y(r1:r2) ones(r2-r1+1,1)];
  end

  % get static and dynamic shims for each slice ................................
  if isempty(magWt), w = ones(lMask,1); else w = magWt(idxMask); end
  switch objNorm % can use different norms for objective function
    case 'L1'
      cvx_begin
        variable shims(size(E,2))
        minimize(norm(w.*(E*shims+axB0(idxMask)), 1))
      cvx_end
    case 'L2'
      %b = E\(-oa(mi));
      shims = lscov(E, -axB0(idxMask), w);
  end
  sShims = shims(1:nsShims); % bs(1):xy, bs(2):x2y2, bs(3):z2
  dShims = zeros(nz, 4); % x, y, z, bo for each slice
  dShims(:,[1 2 4]) = reshape(shims(nsShims+1:end),3,[])';

  % display original and shimmed field maps
  if show
    axB0Shimmed = get_field(xGrid,yGrid,zGrid,sShims,dShims)+axB0;
    imdisp(axB0.*axMask, 'cblab', 'Hz', 'title', 'input field map', ...
      'ucmap', 'jet', 'ulim', [-250 250]);
    imdisp(axB0Shimmed.*axMask, 'cblab', 'Hz', 'title', ...
      'shimmed field (no calibration fix for 2^{nd} order shims)', ...
      'ucmap', 'jet', 'ulim', [-250 250]);
  end

  % adjust z-gradients to minimize (in some sense) z-dephasing .................
  nzSg = size(sgB0,1); % # z pixels in sag map, assume sag map is square
  zRangeSg = sgFov/2*linspace(-1,1,nzSg);
  dz = zRangeSg(2)-zRangeSg(1); % pixel size in z (cm)
  xGridSg = sgXYZ(1)*ones(nzSg);
  % below: not sgXYZ(3)+zRangeSg because fov center moves to isocenter
  [yGridSg,zGridSg] = meshgrid(sgXYZ(2)+zRangeSg, zRangeSg); 
  yVec = vec(yGridSg(1,:)); % y and z axes limits
  zVec = vec(zGridSg(:,1));

  % number of points on each side of slice to do polynomial fit to compute
  % the z shims and update df shims, e.g. npZPfit = 20 means we look +/- 20
  % pixels of the current slice along the z direction.
  npZPfit = round(pfitRangeZ/sgFov*nzSg/2);

  % field shimmed with only second-order gradients
  sgB0SShimmed = sShims(1)*xGridSg.*yGridSg +...
    sShims(2)*(xGridSg.^2-yGridSg.^2) + ...
    sShims(3)*(zGridSg.^2-(xGridSg.^2+yGridSg.^2)/2) + sgB0;

  for i = 1:nz
    % sagittal field map with static shims and dynamic shims for current slice
    sgB0DShimmed = sgMask.*(sgB0SShimmed+dShims(i,1)*xGridSg + ...
      dShims(i,2)*yGridSg + dShims(i,3)*zGridSg + dShims(i,4));

    % find index of axial slice i in the sagittal field map along z direction
    zPos = axXYZ(i,3);
    [~,j] = min(abs(zVec-zPos));

    % here, we build up two vectors zz and ff. zz contains the z positions of
    % all the voxels used for the polynomial fitting for the current slice, and
    % ff contains the dynamic shimmed values (sgB0DShimmed) of these voxels.
    zz = []; ff = [];
    for k = -npZPfit:npZPfit
      % note: the z-location used is the exact obtained from 'ca', but
      % the off-resonance values are from the nearest slice in 'oi'.
      % this shouldn't be a problem if the sagittal resolution is high.
      if j+k>1 && j+k<=nzSg
        oz = sgB0DShimmed(j+k,:); % oz is dummy variable
        oz = oz(sgMask(j+k,:)==1);
        if ~isempty(oz)
          zz = [zz (zPos+dz*k)*ones(1,length(oz))];
          ff = [ff oz];
        end
      end
    end

    % do a polynomial fit over region and update gz and df dynamic shims
    p = polyfit(zz, ff, polyOrd);
    switch polyOrd
      case 1
        dShims(i,3) = dShims(i,3)-p(1);
        dShims(i,4) = dShims(i,4)-p(2);
      case 2
        dShims(i,3) = dShims(i,3)-2*p(1)*zPos-p(2);
        dShims(i,4) = dShims(i,4)+p(1)*zPos^2-p(3);
    end

  end

  if show
    % display original and shimmed field maps
    axB0Shimmed = get_field(xGrid,yGrid,zGrid,sShims,dShims)+axB0;
    imdisp(axB0Shimmed.*axMask, 'cblab', 'Hz', 'title', ...
      'shimmed field (with calibration fix for 2^{nd} order shims)', ...
      'ucmap', 'jet', 'ulim', [-250 250]);
  end

  % account for the static shim imperfections using calibration matrix .........
  coilM = [coilM zeros(9,4)]; % augment matrix 'coilM' to include dynamic shims
  coilM(6,6) = gam;
  coilM(7,7) = gam;
  coilM(8,8) = gam;
  coilM(9,9) = 1;
  
  % shims (static and dynamic) actually needed to generate static shim
  I = coilM\[sShims(1) 0 0 sShims(2) sShims(3) 0 0 0 0]';

  % add these to the per-slice dynamic shims
  gx = I(6)+dShims(:,1)/gam;
  gy = I(7)+dShims(:,2)/gam;
  gz = I(8)+dShims(:,3)/gam;
  dbo = I(9)+dShims(:,4);

  %{
  I   : vector length 5 of second-order amplitudes (A)
  gx  : per-slice Gx (G/cm)
  gy  : per-slice Gy (G/cm)
  gz  : per-slice Gz (G/cm)
  dbo : per-slice freq offset (Hz)
  %}

  % output to shim file ........................................................
  Im = round(1e3*I); % convert to mA
  dShims = [gx gy gz dbo];

  % output info to text file
  if writeFile
    fp = fopen('slshims.txt', 'w');
    fprintf(fp, '(mA)    xy         zy         zx      x2-y2         z2\n');
    fprintf(fp, '%10d %10d %10d %10d %10d\n', ...
      Im(1), Im(2), Im(3), Im(4), Im(5));
    fprintf(fp, '\n');
    fprintf(fp, '(G/cm)     x            y            z      bo (Hz)\n');
    for iz = 1:nz
      fprintf(fp, '%12.6f %12.6f %12.6f %12.6f\n', ...
        gx(iz), gy(iz), gz(iz), dbo(iz));
    end
    fclose(fp);
  end

end

% ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

% extracts the center slices from full dataset
function o = get_centers(O, nH)
  nh = (nH-1)/2;
  o(:,:,1:nh) = O(:,:,2:2:nH);
  o(:,:,nh+1:end) = O(:,:,nH+2:2:end);
end

% ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
  
% returns a grid of positions based on slice center positions
function [xm,ym,zm] = get_xyz(fovx, fovy, nx, ny, c)
  ls = @(x,n) linspace(-x/2,x/2,n);
  [xg yg zm] = meshgrid(ls(fovx,nx), ls(fovx,ny), c(:,3));
  nz = size(c,1);
  xm = zeros(ny, nx, nz);
  ym = zeros(ny, nx, nz);
  % shift x-y center positions slice by slice
  for i = 1:nz
    xm(:,:,i) = c(i,1)+xg(:,:,i);
    ym(:,:,i) = c(i,2)+yg(:,:,i);
  end
end

% ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

% gets the field generated by a set of static and dynamic shims
function o = get_field(x, y, z, bs, bd)
  if nargin == 4, bd = []; end
  nz = size(x,3);
  o = bs(1)*x.*y + bs(2)*(x.^2-y.^2) + bs(3)*(z.^2-(x.^2+y.^2)/2);
  if ~isempty(bd)
    for i = 1:nz
      o(:,:,i) = o(:,:,i)+bd(i,1)*x(:,:,i)+bd(i,2)*y(:,:,i)+...
        bd(i,3)*z(:,:,i)+bd(i,4);
    end
  end
end

% ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

% minimizes ||Ax-b||_1 (L1 minimization)
% code taken from: http://cvxr.com/cvx/examples/html/quickstart.html
function x = l1min(A, b)
  [m,n] = size(A);
  f    = [ zeros(n,1); ones(m,1);  ones(m,1)  ];
  Aeq = sparse(A); Aeq = [Aeq sparse(-eye(m)) sparse(+eye(m)) ];
  %Aeq  = [ A,          -eye(m),    +eye(m)    ];
  lb   = [ -Inf*ones(n,1); zeros(m,1); zeros(m,1) ];
  xzz  = linprog(f,[],[],Aeq,b,lb,[]);
  x = xzz(1:n,:) - xzz(n+1:2*n,:);
end

