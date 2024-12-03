clear
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

T_jacobi = solve_implicit_1d_jacobi(P);
T_direct = solve_implicit_1d_direct(P);
figure(1),plot(P.x,T_direct,'o',P.x,T_jacobi,LineWidth=3);legend('Direct','Jacobi');drawnow

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

P.nx     = 3;                     % Resolution
P.dx     = P.Lx/(P.nx-1);           % Space discretisation step
P.x      = -P.Lx/2:P.dx:P.Lx/2;     % Discretised space
P.S      = P.D/P.dx^2;              % Coefficient in front of T

end

function P = assemble_matrices(P)

nx  = P.nx ;
S   = P.S  ;

% initialisation
A   = zeros(nx,nx);
b   = zeros(nx,1);

for i = 2:nx-1

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
x = zeros(nx, 1);

% Set BC in initial guess. This is necessary, otherwise smoother damping 
% will generate spurious nonzero residuals in Dirichlet nodes
x(1)   = P.TL;
x(end) = P.TR;

P.T_init    = x;

end
 
function T = solve_implicit_1d_jacobi(P)

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

while err >= tol || it >= 1e7
    
    r   = w*inv(diag(D))* (b - A*x);    % Weighted Jacobi iteration
    x   = x + r;                        % Add residual 

    err = norm(x-xold);
    xold = x;
    it  = it + 1;

end

fprintf(" Iterations: %d \n",it)

T = x;

end

function T = solve_implicit_1d_direct(P)

A = P.A;
b = P.b;

res = (A\b);
T   = res';

end