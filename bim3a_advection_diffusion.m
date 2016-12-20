function [A] = bim3a_advection_diffusion(mesh, alpha, beta) 

Nnodes = size(mesh.p,2);
Nelem  = size(mesh.t,2);
Lloc = zeros(8,8,Nelem);
A = sparse(Nnodes,Nnodes);

x = zeros(8,Nelem);
y = zeros(8,Nelem);
z = zeros(8,Nelem);

% coordinate dei punti dei rettangoli
for m = 1:8
    for k = 1:Nelem
    x(m,k) = mesh.p(1, mesh.t(m,k));
    y(m,k) = mesh.p(2, mesh.t(m,k));
    z(m,k) = mesh.p(3, mesh.t(m,k));    
    end
end

% valori hx, hy
hx = x(2,1) - x(1,1);
hy = y(3,1) - y(1,1);
hz = z(5,1) - z(1,1);

%valutazione locale del termine di trasporto
psi12 = zeros(1,Nelem);
psi23 = zeros(1,Nelem);
psi43 = zeros(1,Nelem);
psi14 = zeros(1,Nelem);
psi15 = zeros(1,Nelem);
psi26 = zeros(1,Nelem);
psi37 = zeros(1,Nelem);
psi48 = zeros(1,Nelem);
psi56 = zeros(1,Nelem);
psi67 = zeros(1,Nelem);
psi87 = zeros(1,Nelem);
psi58 = zeros(1,Nelem);
for k = 1:Nelem
        psi12(k) = beta(1,k) * (x(2,k) - x(1,k)) + ...
                   beta(2,k) * (y(2,k) - y(1,k)) + ...
                   beta(3,k) * (z(2,k) - z(1,k));
               
        psi23(k) = beta(1,k) * (x(3,k) - x(2,k)) + ...
                   beta(2,k) * (y(3,k) - y(2,k)) + ...
                   beta(3,k) * (z(3,k) - z(2,k));
        
        psi43(k) = beta(1,k) * (x(3,k) - x(4,k)) + ...
                   beta(2,k) * (y(3,k) - y(4,k)) + ...
                   beta(3,k) * (z(3,k) - z(4,k));
        
        psi14(k) = beta(1,k) * (x(4,k) - x(1,k)) + ...
                   beta(2,k) * (y(4,k) - y(1,k)) + ...
                   beta(3,k) * (z(4,k) - z(1,k));
        
        psi15(k) = beta(1,k) * (x(5,k) - x(1,k)) + ...
                   beta(2,k) * (y(5,k) - y(1,k)) + ...
                   beta(3,k) * (z(5,k) - z(1,k));
        
        psi26(k) = beta(1,k) * (x(6,k) - x(2,k)) + ...
                   beta(2,k) * (y(6,k) - y(2,k)) + ...
                   beta(3,k) * (z(6,k) - z(2,k));
        
        psi37(k) = beta(1,k) * (x(7,k) - x(3,k)) + ...
                   beta(2,k) * (y(7,k) - y(3,k)) + ...
                   beta(3,k) * (z(7,k) - z(3,k));
        
        psi48(k) = beta(1,k) * (x(8,k) - x(4,k)) + ...
                   beta(2,k) * (y(8,k) - y(4,k)) + ...
                   beta(3,k) * (z(8,k) - z(4,k));
        
        psi56(k) = beta(1,k) * (x(6,k) - x(5,k)) + ...
                   beta(2,k) * (y(6,k) - y(5,k)) + ...
                   beta(3,k) * (z(6,k) - z(5,k));
        
        psi67(k) = beta(1,k) * (x(7,k) - x(6,k)) + ...
                   beta(2,k) * (y(7,k) - y(6,k)) + ...
                   beta(3,k) * (z(7,k) - z(6,k));
        
        psi87(k) = beta(1,k) * (x(7,k) - x(8,k)) + ...
                   beta(2,k) * (y(7,k) - y(8,k)) + ...
                   beta(3,k) * (z(7,k) - z(8,k));
        
        psi58(k) = beta(1,k) * (x(8,k) - x(5,k)) + ...
                   beta(2,k) * (y(8,k) - y(5,k)) + ...
                   beta(3,k) * (z(8,k) - z(5,k));
        
end
%valutazione delle bernoulliane
bp12 = zeros(1,Nelem);
bp23 = zeros(1,Nelem);
bp43 = zeros(1,Nelem);
bp14 = zeros(1,Nelem);
bp15 = zeros(1,Nelem);
bp26 = zeros(1,Nelem);
bp37 = zeros(1,Nelem);
bp48 = zeros(1,Nelem);
bp56 = zeros(1,Nelem);
bp67 = zeros(1,Nelem);
bp87 = zeros(1,Nelem);
bp58 = zeros(1,Nelem);

bm12 = zeros(1,Nelem);
bm23 = zeros(1,Nelem);
bm43 = zeros(1,Nelem);
bm14 = zeros(1,Nelem);
bm15 = zeros(1,Nelem);
bm26 = zeros(1,Nelem);
bm37 = zeros(1,Nelem);
bm48 = zeros(1,Nelem);
bm56 = zeros(1,Nelem);
bm67 = zeros(1,Nelem);
bm87 = zeros(1,Nelem);
bm58 = zeros(1,Nelem);

for k = 1:Nelem
    [bp12(k),bm12(k)] = bimu_bernoulli (psi12(k));
    [bp23(k),bm23(k)] = bimu_bernoulli (psi23(k));
    [bp43(k),bm43(k)] = bimu_bernoulli (psi43(k));
    [bp14(k),bm14(k)] = bimu_bernoulli (psi14(k));
    [bp15(k),bm15(k)] = bimu_bernoulli (psi15(k));
    [bp26(k),bm26(k)] = bimu_bernoulli (psi26(k));
    [bp37(k),bm37(k)] = bimu_bernoulli (psi37(k));
    [bp48(k),bm48(k)] = bimu_bernoulli (psi48(k));
    [bp56(k),bm56(k)] = bimu_bernoulli (psi56(k));
    [bp67(k),bm67(k)] = bimu_bernoulli (psi67(k));
    [bp87(k),bm87(k)] = bimu_bernoulli (psi87(k));
    [bp58(k),bm58(k)] = bimu_bernoulli (psi58(k));
    
    bp12(k) = alpha(k) * hy * hz / (4 * hx) * bp12(k);
    bm12(k) = alpha(k) * hy * hz / (4 * hx) * bm12(k);
    bp23(k) = alpha(k) * hx * hz / (4 * hy) * bp23(k);
    bm23(k) = alpha(k) * hx * hz / (4 * hy) * bm23(k);
    bp43(k) = alpha(k) * hy * hz / (4 * hx) * bp43(k);
    bm43(k) = alpha(k) * hy * hz / (4 * hx) * bm43(k);
    bp14(k) = alpha(k) * hx * hz / (4 * hy) * bp14(k);
    bm14(k) = alpha(k) * hx * hz / (4 * hy) * bm14(k);
    
    bp15(k) = alpha(k) * hx * hy / (4 * hz) * bp15(k);
    bm15(k) = alpha(k) * hx * hy / (4 * hz) * bm15(k);
    bp26(k) = alpha(k) * hx * hy / (4 * hz) * bp26(k);
    bm26(k) = alpha(k) * hx * hy / (4 * hz) * bm26(k);
    bp37(k) = alpha(k) * hx * hy / (4 * hz) * bp37(k);
    bm37(k) = alpha(k) * hx * hy / (4 * hz) * bm37(k);
    bp48(k) = alpha(k) * hx * hy / (4 * hz) * bp48(k);
    bm48(k) = alpha(k) * hx * hy / (4 * hz) * bm48(k);
    
    bp56(k) = alpha(k) * hy * hz / (4 * hx) * bp56(k);
    bm56(k) = alpha(k) * hy * hz / (4 * hx) * bm56(k);
    bp67(k) = alpha(k) * hx * hz / (4 * hy) * bp67(k);
    bm67(k) = alpha(k) * hx * hz / (4 * hy) * bm67(k);
    bp87(k) = alpha(k) * hy * hz / (4 * hx) * bp87(k);
    bm87(k) = alpha(k) * hy * hz / (4 * hx) * bm87(k);
    bp58(k) = alpha(k) * hx * hz / (4 * hy) * bp58(k);
    bm58(k) = alpha(k) * hx * hz / (4 * hy) * bm58(k);
end


for k = 1:Nelem
    Lloc(1,1,k) = bm12(k) + bm14(k) + bm15(k);
    Lloc(1,2,k) = -bp12(k);
    Lloc(1,4,k) = -bp14(k);
    Lloc(1,5,k) = -bp15(k);
    
    Lloc(2,1,k) = -bp12(k);
    Lloc(2,2,k) = bm12(k) + bm23(k) + bm26(k); 
    Lloc(2,3,k) = -bp23(k);
    Lloc(2,6,k) = -bp26(k);
    
    Lloc(3,2,k) = -bp23(k);
    Lloc(3,3,k) = bm23(k) + bm43(k) + bm37(k);
    Lloc(3,4,k) = -bp43(k);
    Lloc(3,7,k) = -bp37(k);
    
    Lloc(4,1,k) = -bp14(k);
    Lloc(4,3,k) = -bp43(k);
    Lloc(4,4,k) = bm43(k) + bm14(k) + bm48(k);
    Lloc(4,8,k) = -bp48(k);
    
    Lloc(5,1,k) = -bp15(k);
    Lloc(5,5,k) = bm15(k) + bm56(k) + bm58(k);
    Lloc(5,6,k) = -bp56(k);
    Lloc(5,8,k) = -bp58(k);
    
    Lloc(6,2,k) = -bp26(k);
    Lloc(6,5,k) = -bp56(k);
    Lloc(6,6,k) = bm56(k) + bm26(k) + bm67(k);
    Lloc(6,7,k) = -bp67(k);
    
    Lloc(7,3,k) = -bp37(k);
    Lloc(7,6,k) = -bp67(k);
    Lloc(7,7,k) = bm37(k) + bm67(k) + bm87(k);
    Lloc(7,8,k) = -bp87(k);
    
    Lloc(8,4,k) = -bp48(k);
    Lloc(8,5,k) = -bp58(k);
    Lloc(8,7,k) = -bp87(k);
    Lloc(8,8,k) = bm48(k) + bm58(k) + bm87(k);
    
    for local_i = 1:8
        for local_j = 1:8
            A(mesh.t(local_i,k), mesh.t(local_j,k)) = ...
                A(mesh.t(local_i,k), mesh.t(local_j,k)) + Lloc(local_i,local_j,k); 
        end
    end
end

end