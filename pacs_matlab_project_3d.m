%clc
%L2 = [];
%H1 = [];
%for NH = 1:10
%HXX = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001,0.0005,0.0002,0.0001];
hx = .1;%HXX(NH);
hy = 1/3;
hz = 1/3;

mesh = msh3m_quadmesh(0:hx:1, 0:hy:1, 0:hz:1, 1, 1:6);


x    = mesh.p(1,:)';

Dnodes  = msh3m_nodes_on_faces(mesh,[1,2]);
Nnodes = size(mesh.p, 2);
Nelements = size(mesh.t,2);
Varnodes = setdiff(1:Nnodes, Dnodes);

alpha = 0.1*ones(Nelements,1);
%eta    = .1*ones(Nnodes,1);
beta = [ones(1,Nelements)/0.1;
        0*ones(1,Nelements)/0.1;
        0*ones(1,Nelements)/0.1];
%gamma  = ones(Nnodes,1);
f      = bim3a_rhs(mesh, ones (Nelements, 1), ones (Nnodes, 1));

S = bim3a_advection_diffusion(mesh,alpha,beta);
u = zeros(Nnodes,1);
uex = x - (exp(10*x)-1)/(exp(10)-1);
u(Varnodes) = S(Varnodes,Varnodes) \ f(Varnodes);
%assert(max(abs(u-uex)) < 1e-7)

Nx = 1/hx; Ny = 1/hy; Nz = 1/hz;
X = zeros(Ny+1,Nx+1);
for i = 1:Nx+1
    X(:,i) = u( (i-1)*(Ny+1)+1 : i*(Ny+1));
end
surf(0:hx:1, 0:hy:1, X )

%u_ex = @(x,y) x - (exp(10*x)-1)/(exp(10)-1) + 0*y;

%dxuex = @(x,y) 1 - ( 10*exp(10*x) )/(exp(10)-1) +0*y;
%dyuex = @(x,y) 0*x+0*y;
%[L2error, H1error] = rectFE_error(mesh,u_ex,dxuex,dyuex,u);
%[L2] = [L2,L2error];
%[H1] = [H1,H1error];
%end
