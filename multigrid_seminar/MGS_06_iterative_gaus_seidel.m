% 101x +  12y - 13z  = 14
%  21x + 201y + 23z  = 24
% -31x +  32y + 301z = 34
%  Ax = b 
%  x = A^(-1)*b

% Assemble matrices and vectors 0.142996780344579	0.0909589332729228	0.118014001082885
% Matrix A must be diagonally dominant! 
A = [101 12  -13  ;...
      21 201 23   ;...
     -31 32  301 ];

b = [ 14 ;
      24 ;
      34 ];

x = [ 1e3  ;
      -2e3  ;
      5e6 ];

% initalise values
err  = 1e10;    % high initial error to update
tol  = 1e-9;    % preferred accuracy for solution
it   = 0;       % iterations counter
xold = x;
D    = diag(A);           % get the diagonal of matrix A
L    = tril(A); % get lower triangular part of matrix A
U    = triu(A) - diag(D); % get upper triangular part of matrix A

while err >= tol || it >= 1e7

    x   = inv(L)* (b - (U)*x);
    err = norm(x-xold);
    xold = x;
    it  = it + 1;
    
end

% Get results from direct
xan  = A\b;
fprintf("Analytical solution is\n")
fprintf( " %g \n",xan)
fprintf("yours is \n" )
fprintf("%g \n ", x )
fprintf(" Iterations: %d \n",it)


