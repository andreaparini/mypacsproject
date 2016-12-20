function msh = msh3m_quadmesh(x,y,z,region,sides)

  nx = length (x);
  ny = length (y);
  nz = length (z);

  [XX, YY, ZZ] = meshgrid (x, y, z);
  msh.p = [XX(:), YY(:), ZZ(:)]';  
  
  iiv (ny,nx,nz)=0;
  iiv(:)=1:nx*ny*nz;
  iiv(end,:,:)=[];
  iiv(:,end,:)=[];
  iiv(:,:,end)=[];
  iiv=iiv(:)';
  N1 = nx*ny;
  
  msh.t = [iiv;    iiv+ny;    iiv+ny+1;    iiv+1;
           iiv+N1; iiv+ny+N1; iiv+ny+1+N1; iiv+1+N1];
  
  
  % generate boundary face list
  
  % left
  T = msh.t;
  T(:) = msh.p(1, msh.t)' == x(1);
  [~, order] = sort (T, 1);
  ii = (find(sum(T,1)==4));
  order(1:4,:) = [];
  for jj=1:length (ii)
    e1(:,jj) = msh.t(order(:,ii(jj)),ii(jj));
  end
  e1(10,:) = sides(1);
  
    % right
  T(:) = msh.p(1,msh.t)' == x(end);
  [~, order] = sort (T, 1);
  ii = (find (sum (T, 1) == 4));
  order(1:4,:) = [];
  for jj=1:length (ii)
    e2(:,jj) = msh.t(order(:,ii(jj)),ii(jj));
  end
  e2(10,:) = sides(2);

  % front
  T(:) = msh.p(2,msh.t)' == y(1);
  [~, order] = sort (T, 1);
  ii = (find (sum (T, 1) == 4));
  order(1:4,:) = [];
  for jj=1:length (ii)
    e3(:,jj) = msh.t(order(:,ii(jj)),ii(jj));
  end
  e3(10,:) = sides(3);

  % back
  T(:) = msh.p(2,msh.t)' == y(end);
  [~,order] = sort (T,1);
  ii = (find (sum (T,1) == 4));
  order(1:4,:) = [];
  for jj=1:length (ii)
    e4(:,jj) = msh.t(order(:,ii(jj)),ii(jj));
  end
  e4(10,:) = sides(4);
  
  % bottom
  T       = msh.t;
  T(:)    = msh.p(3,msh.t)'==z(1);
  [~,order] = sort(T,1);
  ii      = (find (sum (T,1)==4));
  order(1:4,:) = [];
  for jj=1:length(ii)
    e5(:,jj)      = msh.t(order(:,ii(jj)),ii(jj));
  end
  e5(10,:) = sides(5);
  
  % top
  T       = msh.t;
  T(:)    = msh.p(3,msh.t)'==z(end);
  [~,order] = sort(T,1);
  ii      = (find (sum (T,1) == 4));
  order(1:4,:) = [];
  for jj=1:length(ii)
    e6(:,jj)      = msh.t(order(:,ii(jj)),ii(jj));
  end
  e6(10,:) = sides(6);

  msh.e       = [e1,e2,e3,e4,e5,e6];
  msh.e (11,:) = region;
  msh.t(9,:) = region;
end