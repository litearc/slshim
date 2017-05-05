function M = coilm(o, m, A, zp, varargin)
  %
  %  calculates the coil calibration matrix M [9 x 5], where M(i,j) is the
  %  amplitude of field i produced by coil j. we model each second-order coil j
  %  as a sum of the fields i: xy zy zx x2-y2 z2-(x2+y2)/2 x y z bo. the
  %  amplitudes have units of Hz/(cm2*A), Hz/(cm*A), or Hz/A.
  %
  %  function M = coilm(o, m, A, zp)
  %
  %  inputs ....................................................................
  %  o                cell array of field maps {baseline, xy, zy, zx, x2y2, z2}.
  %  m                mask. [x y z] (same size as field maps)
  %  A                # amps input to each shim coil (above baseline). (float)
  %  zp               slice z-positions (assume x and y at center). (vector)
  %
  %  outputs ...................................................................
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

  [dte,fov,tell] = setopts(varargin, {'dte',1,'fov',22,'tell',0});

  % subtract baseline
  for i = 2:6
    o{i} = m.*(o{i}-o{1});
  end

  % get position vectors
  v1 = linspace(-fov/2, fov/2, size(m,1));
  [x, y, z] = meshgrid(v1, v1, zp);
  ii = find(m==1); x = x(ii); y = y(ii); z = z(ii);

  % field-spatial dependence functions: bo(x,y,z) = gam * G{i} * f{i}(x,y,z)
  fxy = x.*y;
  fzy = z.*y;
  fzx = z.*x;
  fx2y2 = x.^2-y.^2;
  fz2 = z.^2-.5*(x.^2+y.^2);

  % do inversion - units are: Hz/cm2/A for 2nd order, Hz/cm/A for linear, Hz/A for offset
  M = zeros(9, 5);
  for i = 2:6
    ga = [fxy fzy fzx fx2y2 fz2 x y z ones(length(x),1)] \ o{i}(ii);
    M(:,i-1) = ga/A;
  end

  % calculate percent variance explained for each field
  if tell
    [x, y, z] = meshgrid(v1, v1, zp);
    fxy = x.*y;
    fzy = z.*y;
    fzx = z.*x;
    fx2y2 = x.^2-y.^2;
    fz2 = z.^2-.5*(x.^2+y.^2);

    mf = {}; % modeled field
    ii = find(m);
    names = {'xy','zy','zx','x2y2','z2'};
    vara = @(x) var(x(:));
    for i = 1:5
      mf{i} = M(1,i)*fxy + M(2,i)*fzy + M(3,i)*fzx + M(4,i)*fx2y2 + M(5,i)*fz2 + ...
        M(6,i)*x + M(7,i)*y + M(8,i)*z + M(9,i);
      f1 = mf{i}*A;
      f2 = o{i+1};
      imdisp(f1.*m-f2, 'ucmap', 'jet', 'cblab', 'Hz');
      ef = 1-vara(f2(ii)-f1(ii))/vara(f2(ii));
      fprintf('exp var %s: %f\n', names{i}, ef);
    end
  end

end
