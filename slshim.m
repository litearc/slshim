function [I bd] = slshim(oa, ma, ca, fova, os, ms, ss, fovs, M, varargin)
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
  %  function [I bd] = slshim(oa, ma, ca, sa, fova, os, ms, cs, ss, fovs, M)
  %
  %  inputs ....................................................................
  %  oa               axial B0 map. [x y z] (Hz)
  %  ma               mask for `oa'. [x y z] (binary)
  %  ca               center position (x,y,z) for each slice in `oa'.
  %                   [slices (x,y,z)] (cm)
  %                   the z-locations are centered around 0 due to the table
  %                   shifting
  %  fova             field of view for `oa'. (cm)
  %  os               sagittal B0 map, 1 slice. [z x] (Hz)
  %  ms               mask for `os'. [z x]
  %  ss               position (x,y,z) of slice in `os' (cm)
  %  fovs             field of view for 'os' (cm)
  %  M                matrix where M(i,j) is amount of field i that coil j
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
  %  out              output shim file? (default = 0)
  %  gam              gamma (Hz/G) (default = 4258)
  %  obj_norm         norm to use for objective function
  %                   options: 'L1', 'L2' (default)
  %  magw             mag image used as weights in optimization
  %  pfdz             thickness in z (cm) used to do polynomial fit for
  %                   z-gradient calculation (default = 8)
  %  ord              polynomial order for z-gradient calculation (default = 2)
  %  sag              if 1, display sagittal mag image of brain with slices
  %  scm              spinal-cord mask (default = none)
  %  show             display various images (default = 0)
  %
  %  outputs ...................................................................
  %  I                second-order amplitudes: [xy zy zx x2-y2 z2-.5(x2+y2)]
  %  bd               linear + offset for each slice: [x y z o] (nslices x 3)
  %
 
  [out, gam, obj_norm, magw, pfdz, ord, sag, scm, show] = setopts(varargin, {...
    'out', 0, 'gam', 4258, 'obj_norm', 'L2', 'magw', [], 'pfdz', 8, ...
    'ord', 2, 'sag', [], 'scm', [], 'show', 0});

  % ............................................................................
  [ny, nx, nz] = size(oa);
  ns = 3; % # static shims (excluding xz and yz)
  nd = 3; % # dynamic shims (excluding z, which is handled later)
  oa = oa.*ma; % apply masks
  os = os.*ms;

  [xm,ym,zm] = get_xyz(fova, fova, nx, ny, ca);
 
  % for optimization, we only care about the field inside mask
  mi = find(ma); l = length(mi); x = xm(mi); y = ym(mi); z = zm(mi);

  % 'in' contains the index of the first pixel inside the mask for each slice
  in = zeros(1,nz); for i = 1:nz, in(i) = find(z==ca(i,3),1); end

  % position encoding matrices for static and dynamic shims ....................
  Es = [x.*y, x.^2-y.^2, z.^2-(x.^2+y.^2)/2];
  E = [Es, zeros(l,nd*nz)];
  % for each slice, fill in one "block" section in the matrix
  % c1,c2 are the beginning/end columns
  % r1,r2 are the beginning/end rows
  for i = 1:nz
    c1 = ns+1+nd*(i-1);
    c2 = ns+nd*i;
    r1 = in(i);
    if i == nz, r2 = l; else r2 = in(i+1)-1; end
    E(r1:r2, c1:c2) = [x(r1:r2) y(r1:r2) ones(r2-r1+1,1)];
  end

  % get static and dynamic shims for each slice ................................
  if isempty(magw), w = ones(l,1); else w = magw(mi); end
  switch obj_norm % can use different norms for objective function
    case 'L1'
      cvx_begin
        variable b(size(E,2))
        minimize(norm(w.*(E*b+oa(mi)), 1))
      cvx_end
    case 'L2'
      %b = E\(-oa(mi));
      b = lscov(E, -oa(mi), w);
  end
  bs = b(1:ns); % bs(1):xy, bs(2):x2y2, bs(3):z2
  bd = zeros(nz, 4); % x, y, z, bo for each slice
  bd(:,[1 2 4]) = reshape(b(ns+1:end),3,[])';

  % display original and shimmed field maps
  if show
    oo = get_field(xm,ym,zm,bs,bd)+oa;
    imdisp(oa.*ma, 'cblab', 'Hz', 'title', 'input field map', 'ucmap', 'jet', ...
      'ulim', [-250 250]);
    imdisp(oo.*ma, 'cblab', 'Hz', 'title', ...
      'shimmed field (no calibration fix for 2^{nd} order shims)', ...
      'ucmap', 'jet', 'ulim', [-250 250]);
  end

  % adjust z-gradients to minimize (in some sense) z-dephasing .................
  np = size(os,1); % assume square
  sz = fovs/2*linspace(-1,1,np);
  dz = sz(2)-sz(1); % pixel size in z (cm)
  xs = ss(1)*ones(np);
  [ys,zs] = meshgrid(ss(2)+sz, sz); % not ss(3)+sz because fov center moves to isocenter
  yv = vec(ys(1,:)); % y and z axes limits
  zv = vec(zs(:,1));
  npzc = round(pfdz/fovs*np/2);

  % field shimmed with only second-order gradients
  oo = bs(1)*xs.*ys+bs(2)*(xs.^2-ys.^2)+bs(3)*(zs.^2-(xs.^2+ys.^2)/2)+os;
  oos = zeros(np); % sagittal image where slices collected are shimmed values

  % display axial and sagittal field maps
  if show
    minv = min([oa(:); os(:)]);
    maxv = max([oa(:); os(:)]);
    figure; hold on;
    imagesc(yv, zv, os);
    axis square; colormap jet; cbar('label','Hz');
    title('input field map (sagittal)');
    for i = 1:nz
      plot([yv(1) yv(end)], ca(i,3)*[1 1], 'k');
    end
    xlim([yv(1) yv(end)]);
    ylim([zv(1) zv(end)]);
    set(gca, 'ydir', 'reverse');
    caxis([-250 250]);
    % caxis([minv maxv]);
  end

  for i = 1:nz
    % sagittal field map for current slice
    oi = ms.*(oo+bd(i,1)*xs+bd(i,2)*ys+bd(i,3)*zs+bd(i,4));

    % find closest pixel in z dimension for current slice
    [~,j] = min(abs(zv-ca(i,3)));

    zp = ca(i,3);
    xx = []; yy = [];
    for k = -npzc:npzc
      % note: the z-location used is the exact obtained from 'ca', but
      % the off-resonance values are from the nearest slice in 'oi'.
      % this shouldn't be a problem if the sagittal resolution is high.
      if j+k>1 && j+k<=np
        oz = oi(j+k,:);
        oz = oz(ms(j+k,:)==1);
        if ~isempty(oz)
          xx = [xx (zp+dz*k)*ones(1,length(oz))];
          yy = [yy oz];
        end
      end
    end

    % do a polynomial fit over region and update gz and df dynamic shims
    p = polyfit(xx, yy, ord);
    xvv = linspace(xx(1),xx(end),256);
    switch ord
      case 1
        bd(i,3) = bd(i,3)-p(1);
        bd(i,4) = bd(i,4)-p(2);
        oos(j,:) = oi(j,:)-p(1)*zp-p(2);
      case 2
        bd(i,3) = bd(i,3)-2*p(1)*zp-p(2);
        bd(i,4) = bd(i,4)+p(1)*zp^2-p(3);
        oos(j,:) = oi(j,:)-p(1)*zp^2-p(2)*zp-p(3);
    end

    % for thesis presentation
    if i == 22 && show
      figure; fig(8,6); hold on;
      xlim(lims(xx));
      xlabel('z (cm)'); ylabel('\Deltaf (Hz)');
      ax = plot(xx, yy, 'b.'); set(ax, 'Color', [0 .4 .8]);
      ax = plot(xvv, p(1)*xvv.^2+p(2)*xvv+p(3), 'linewidth', 2);
      set(ax, 'Color', [.8 .2 .2]);
      ax = plot(zp, p(1)*zp^2+p(2)*zp+p(3), 'k.');
      set(ax, 'markersize', 20);
      ax = plot(xvv, (2*p(1)*zp+p(2))*xvv - (2*p(1)*zp+p(2))*zp + (p(1)*zp.^2+p(2)*zp+p(3)), 'k', 'linewidth', 2);
      set(gca, 'position', [.2 .2 .6 .7]);
      title('z and df shim calculation for slice 22');
    end

  end

  if show
    % display original and shimmed field maps
    oo = get_field(xm,ym,zm,bs,bd)+oa;
    imdisp(oo.*ma, 'cblab', 'Hz', 'title', ...
      'shimmed field (with calibration fix for 2^{nd} order shims)', ...
      'ucmap', 'jet', 'ulim', [-250 250]);
  end

  % account for the static shim imperfections using calibration matrix .........
  M = [M zeros(9,4)]; % augment matrix 'M' to include dynamic shims
  M(6,6) = gam;
  M(7,7) = gam;
  M(8,8) = gam;
  M(9,9) = 1;
  
  % shims (static and dynamic) actually needed to generate static shim
  I = M\[bs(1) 0 0 bs(2) bs(3) 0 0 0 0]';

  % add these to the per-slice dynamic shims
  gx = I(6)+bd(:,1)/gam;
  gy = I(7)+bd(:,2)/gam;
  gz = I(8)+bd(:,3)/gam;
  dbo = I(9)+bd(:,4);

  %{
  I   : vector length 5 of second-order amplitudes (A)
  gx  : per-slice Gx (G/cm)
  gy  : per-slice Gy (G/cm)
  gz  : per-slice Gz (G/cm)
  dbo : per-slice freq offset (Hz)
  %}

  % output to shim file ........................................................
  Im = round(1e3*I); % convert to mA
  bd = [gx gy gz dbo];

  % output info to text file
  if out
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

