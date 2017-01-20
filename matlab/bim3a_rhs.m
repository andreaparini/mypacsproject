function rhs = bim3a_rhs(mesh,f,g)

Nnodes = size(mesh.p,2);
Nelem  = size(mesh.t,2);
rhs = sparse(Nnodes,1);
x = zeros(8,Nelem);
y = zeros(8,Nelem);
z = zeros(8,Nelem);
F = zeros(8,Nelem);
% coordinate dei punti dei rettangoli e valutazione locale di f
for m = 1:8
    for k = 1:Nelem
    x(m,k) = mesh.p(1, mesh.t(m,k));
    y(m,k) = mesh.p(2, mesh.t(m,k));
    z(m,k) = mesh.p(3, mesh.t(m,k));
    F(m,k) = f(k) * g(mesh.t(m,k));
    end
end

% valori hx, hy
hx = x(2,1) - x(1,1);
hy = y(3,1) - y(1,1);
hz = z(5,1) - z(1,1);

rhs_loc = zeros(8,Nelem);
for k = 1:Nelem
    rhs_loc(1,k) = hx * hy * hz / 8 * F(1,k);
    rhs_loc(2,k) = hx * hy * hz / 8 * F(2,k);
    rhs_loc(3,k) = hx * hy * hz / 8 * F(3,k);
    rhs_loc(4,k) = hx * hy * hz / 8 * F(4,k);
    rhs_loc(5,k) = hx * hy * hz / 8 * F(5,k);
    rhs_loc(6,k) = hx * hy * hz / 8 * F(6,k);
    rhs_loc(7,k) = hx * hy * hz / 8 * F(7,k);
    rhs_loc(8,k) = hx * hy * hz / 8 * F(8,k);
    for local_j = 1:8
            rhs(mesh.t(local_j,k)) = ...
                rhs(mesh.t(local_j,k)) + rhs_loc(local_j,k); 
    end
end

end