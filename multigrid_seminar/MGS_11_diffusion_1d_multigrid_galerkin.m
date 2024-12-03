clear,clc, figure(1),figure(2),clf(1)
% 1D Thermal diffusion
% D*d2T/dx2=0
% D*(T(i-1)-2*T(i)+T(i+1))/dx^2=0
% S = D/dx^2
% S - 2S + S = 0
% A = [1  0  0
%      S -2S S
%      0  0  1 ]

P = [];
P = initialise_physic_const(P); % Density, conductivity, length etc
P = intialise_numeric_const(P); % Resolution, discretised space
P = assemble_matrices(P);
P = assemble_multigrid(P);

P.T_MG     = solve_implicit_1d_multigrid(P);
P.T_direct = solve_implicit_1d_direct(P);

figure(1),plot(P.x,P.T_direct,'o',P.x,P.T_MG,LineWidth=3);legend('Direct','Multigrid');drawnow;

T_jacobi   = solve_implicit_1d_jacobi(P);

figure(1),plot(P.x,P.T_direct,'o',P.x,P.T_MG,P.x,T_jacobi,'|',LineWidth=3);legend('Direct','Multigrid','Jacobi');drawnow;


function P = initialise_physic_const(P)

P.Lx     = 10;              % Length of the domain
P.TL     = 1;               % Temperature at left boundary
P.TR     = 0;               % Temperature at right boundary
P.k      = 1;               % Conductivity
P.rho    = 1;               % density
P.Cp     = 1;               % Specific heat capacity
P.D      = P.k/P.rho/P.Cp;  % Diffusivity

end

function P = intialise_numeric_const(P)

P.nc           = 2^3;                    % Number of cells
P.n            = P.nc + 1 ;              % Number of nodes
P.dx           = P.Lx/(P.n-1);           % Space discretisation step
P.x            = -P.Lx/2:P.dx:P.Lx/2;    % Discretised space
P.S            = P.D/P.dx^2;             % Coefficient in front of T
P.nlevels      = 3;                     % Number of multigrid layers
P.nsweeps      = 10;                     % Number of multigrid iterations
P.maxit_smooth = 5;                      % Number of iterations during smoothing

end

function P = assemble_matrices(P)

n   = P.n ;
S   = P.S  ;

% initialisation
A   = zeros(n,n);
b   = zeros(n,1);

for i = 2:n-1

    A(i,i-1)  =    -S;
    A(i,i  )  =   2*S;
    A(i,i+1)  =    -S;

end

% LEFT Boundary
A(1,1)     = 1;
b(1)       = P.TL;

% Right Boundary
A(end,end) = 1;
b(end)     = P.TR;

P.A = A;
P.b = b;

% initial guess
x = zeros(n, 1);

% Set BC in initial guess. This is necessary, otherwise smoother damping 
% will generate spurious nonzero residuals in Dirichlet nodes
x(1)   = P.TL;
x(end) = P.TR;

P.T_init = x;

end

function P = assemble_multigrid(P)

    inc = P.nc;

    for i = 1:P.nlevels
        
        % get restriction matrix for all levels but last
        if i ~= P.nlevels
            
            P.levels(i).R = getRestrict(inc);
             
            inc = inc/2;

        end
    
        % get linear operator
        if i == 1
            
            % top level
            P.levels(i).A = P.A;

        else
            
            % Galerkin coarsening
            P.levels(i).A = P.levels(i-1).R * P.levels(i-1).A * P.levels(i-1).R';

        end

        P.levels(i).invDiag  = 1./(full(diag(P.levels(i).A)));
        
    end

end

function [R] = getRestrict(nc)

% uniform restriction by factor of two

if mod(nc, 2) ~= 0
    error("Odd number of cells specified")
end
    
% number of cells in restricted grid
ncr = nc/2;

% number of nodes
n  = nc  + 1;
nr = ncr + 1;

R = sparse(nr, n);

for i = 2:nr-1
    
    R(i, 2*i  ) =  0.5;
    R(i, 2*i-1) =  1.0;
    R(i, 2*i-2) =  0.5;

end

% BC at left side
R(1, 1) = 1.0;

% BC at right side
R(nr, n) = 1.0;

end

function [x] = vcycle(P,T)

levels  = P.levels;
nlevels = length(levels);
maxit   = P.maxit_smooth;
b       = P.b;
x       = T;

% pre-smooth solution
x = smoother(levels(1).A, b, x, levels(1).invDiag, maxit);

% get top level residual
r = b - levels(1).A*x;

% restrict residual to next level
levels(2).r = levels(1).R*r;

% perform restriction step
for i = 2:nlevels-1
    
    % apply pre-smoothing
    levels(i).e = smoother(levels(i).A, levels(i).r, zeros(size(levels(i).r)),levels(i).invDiag, maxit);
       
    % update residual
    levels(i).r = levels(i).r - levels(i).A*levels(i).e;
          
    % restrict residual to next level
    levels(i+1).r = levels(i).R*levels(i).r;
    
end

% coarse grid correction
levels(nlevels).e = levels(nlevels).A\levels(nlevels).r;

% perform prolongation step
for i = nlevels-1:-1:2
    
    % add correction from previous level
    levels(i).e = levels(i).e + levels(i).R'*levels(i+1).e; 
    
    % apply post-smoothing
    levels(i).e = smoother(levels(i).A, levels(i).r, levels(i).e,levels(i).invDiag, maxit);
    
end

% update solution
x = x + levels(1).R'*levels(2).e;

% post-smooth solution
x = smoother(levels(1).A, b, x,levels(1).invDiag, maxit);

end

function [x] = smoother(A, b, x,invDiag, maxit)

% Richardson iteration with Jacobi preconditioner

w    = 2/3;

for i = 1:maxit

    r   = w.* (b - A*x).*invDiag;       % Weighted Jacobi iteration
    x   = x + r;                        % Add residual

end

end

function T = solve_implicit_1d_multigrid(P)

t1 = tic;
T  = P.T_init;

for i = 1:P.nsweeps

    T = vcycle(P,T);

end

t2 = toc(t1);
fprintf('Multigrid solver done in %f seconds \n', t2)

end

function T = solve_implicit_1d_jacobi(P)

t1 = tic;

A = P.A;
b = P.b;
x = P.T_init;

% initalise values
err  = 1e10;    % high initial error to update
tol  = 1e-6;    % preferred accuracy for solution
it   = 0;       % iterations counter
xold = x;
D    = diag(A);           % get the diagonal of matrix A
w    = 2/3;
maxit = 1e5;

while err >= tol || it >= maxit
    
    r   = w* (b - A*x).* P.levels(1).invDiag;    % Weighted Jacobi iteration
    x   = x + r;                        % Add residual 

    err = norm(x-xold);
    xold = x;
    it  = it + 1;
    if mod(it,round(P.n^(1/2))) == 0

%     figure(1),plot(P.x,P.T_direct,'o',P.x,P.T_MG,P.x,x,'|',LineWidth=3);legend('Direct','Multigrid','Jacobi');drawnow;

    end
end

T = x;

t2 = toc(t1);
if it>=maxit; fprintf(" Maximum iterations reached: %d \n",it);end;
fprintf('Jacobi    solver done in %f seconds \n', t2)

end

function T = solve_implicit_1d_direct(P)

t1 = tic;

A = P.A;
b = P.b;

res = (A\b);
T   = res';

t2  = toc(t1);
fprintf('Direct    solver done in %f seconds \n', t2)

end