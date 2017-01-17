function rhs = bim2a_rhs(mesh,f,g)

Nnodes = size(mesh.p,2);
Nelem  = size(mesh.t,2);
rhs = sparse(Nnodes,1);
x = zeros(4,Nelem);
y = zeros(4,Nelem);
F = zeros(4,Nelem);
% coordinate dei punti dei rettangoli e valutazione locale di f
for m = 1:4
    for k = 1:Nelem
    x(m,k) = mesh.p(1, mesh.t(m,k));
    y(m,k) = mesh.p(2, mesh.t(m,k));
    F(m,k) = f(k) * g(mesh.t(m,k));
    end
end

% valori hx, hy
hx = x(2,1) - x(1,1);
hy = y(3,1) - y(1,1);

rhs_loc = zeros(4,Nelem);
for k = 1:Nelem
    rhs_loc(1,k) = hx * hy / 4 * F(1,k);
    rhs_loc(2,k) = hx * hy / 4 * F(2,k);
    rhs_loc(3,k) = hx * hy / 4 * F(3,k);
    rhs_loc(4,k) = hx * hy / 4 * F(4,k);
    for local_j = 1:4
            rhs(mesh.t(local_j,k)) = ...
                rhs(mesh.t(local_j,k)) + rhs_loc(local_j,k); 
    end
end

end