function [A] = bim2a_advection_diffusion(mesh, alpha, beta) 

Nnodes = size(mesh.p,2);
Nelem  = size(mesh.t,2);
Lloc = zeros(4,4,Nelem);
A = sparse(Nnodes,Nnodes);

x = zeros(4,Nelem);
y = zeros(4,Nelem);

% coordinate dei punti dei rettangoli
for m = 1:4
    for k = 1:Nelem
    x(m,k) = mesh.p(1, mesh.t(m,k));
    y(m,k) = mesh.p(2, mesh.t(m,k));
    end
end

% valori hx, hy
hx = x(2,1) - x(1,1);
hy = y(3,1) - y(1,1);

%valutazione locale del termine di trasporto
for k = 1:Nelem
        psi12(k) = beta(1,k) * (x(2,k) - x(1,k)) + ...
                   beta(2,k) * (y(2,k) - y(1,k));
        psi23(k) = beta(1,k) * (x(3,k) - x(2,k)) + ...
                   beta(2,k) * (y(3,k) - y(2,k));
        psi34(k) = beta(1,k) * (x(4,k) - x(3,k)) + ...
                   beta(2,k) * (y(4,k) - y(3,k));
        psi41(k) = beta(1,k) * (x(1,k) - x(4,k)) + ...
                   beta(2,k) * (y(1,k) - y(4,k));
end
%valutazione delle bernoulliane
bp12 = zeros(1,Nelem);
bp23 = zeros(1,Nelem);
bp34 = zeros(1,Nelem);
bp41 = zeros(1,Nelem);
bm12 = zeros(1,Nelem);
bm23 = zeros(1,Nelem);
bm34 = zeros(1,Nelem);
bm41 = zeros(1,Nelem);
for k = 1:Nelem
    [bp12(k),bm12(k)] = bimu_bernoulli (psi12(k));
    [bp23(k),bm23(k)] = bimu_bernoulli (psi23(k));
    [bp34(k),bm34(k)] = bimu_bernoulli (psi34(k));
    [bp41(k),bm41(k)] = bimu_bernoulli (psi41(k));
    
    bp12(k) = alpha(k) * hy / (2 * hx) * bp12(k);
    bm12(k) = alpha(k) * hy / (2 * hx) * bm12(k);
    bp23(k) = alpha(k) * hx / (2 * hy) * bp23(k);
    bm23(k) = alpha(k) * hx / (2 * hy) * bm23(k);
    bp34(k) = alpha(k) * hy / (2 * hx) * bp34(k);
    bm34(k) = alpha(k) * hy / (2 * hx) * bm34(k);
    bp41(k) = alpha(k) * hx / (2 * hy) * bp41(k);
    bm41(k) = alpha(k) * hx / (2 * hy) * bm41(k);
end


for k = 1:Nelem
    Lloc(1,1,k) = bm12(k) + bp41(k);
    Lloc(1,2,k) = -bp12(k);
    Lloc(1,4,k) = -bm41(k);
    
    Lloc(2,1,k) = -bm12(k);
    Lloc(2,2,k) = bp12(k) + bm23(k);
    Lloc(2,3,k) = -bp23(k);
   
    Lloc(3,2,k) = -bm23(k);
    Lloc(3,3,k) = bp23(k) + bm34(k);
    Lloc(3,4,k) = -bp34(k);
    
    Lloc(4,1,k) = -bp41(k);
    Lloc(4,3,k) = -bm34(k);
    Lloc(4,4,k) = bp34(k) + bm41(k);
    for local_i = 1:4
        for local_j = 1:4
            A(mesh.t(local_i,k), mesh.t(local_j,k)) = ...
                A(mesh.t(local_i,k), mesh.t(local_j,k)) + Lloc(local_i,local_j,k); 
        end
    end
end

end