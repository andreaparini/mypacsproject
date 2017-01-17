function msh = msh2m_quadmesh (x, y, region, sides)
  
  nx = numel (x);
  ny = numel (y);
  
  [XX,YY] = meshgrid (x, y);
  msh.p = [XX(:), YY(:)].';

  iiv(ny,nx) = 0;
  iiv(:) = 1 : nx*ny;
  iiv(end,:) = [];
  iiv(:,end) = [];
  iiv = iiv(:).';

  msh.t = [iiv; iiv+ny; iiv+ny+1; iiv+1];
  msh.t(5,:) = region;

  l1 = 1 + ny * ([1:nx] - 1);
  l4 = 1:ny;
  l2 = ny * (nx-1) + 1:nx*ny;
  l3 = ny + l1 - 1;

  [s1, s2, s3, s4] = deal (sides(1), sides(2), sides(3), sides(4));
  msh.e = [l1([1:end-1])      l2([1:end-1])      l3([1:end-1])      l4([1:end-1])
           l1([2:end])        l2([2:end])        l3([2:end])        l4([2:end])
           [l1([1:end-1])     l2([1:end-1])      l3([1:end-1])      l4([1:end-1])]*0
           [l1([1:end-1])     l2([1:end-1])      l3([1:end-1])      l4([1:end-1])]*0
           l1([1:end-1])*0+s1 l2([1:end-1])*0+s2 l3([1:end-1])*0+s3 l4([1:end-1])*0+s4
           [l1([1:end-1])     l2([1:end-1])      l3([1:end-1])      l4([1:end-1])]*0
           [l1([1:end-1])     l2([1:end-1])      l3([1:end-1])      l4([1:end-1])]*0+region];
  
end %function
