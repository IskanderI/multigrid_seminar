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

x = [ 0  ;
      0  ;
      0 ];

% Initialise values
err  = 1e10;    % high initial error to update
tol  = 1e-7;    % preferred accuracy for solution
it   = 0;       % iterations counter
xold = x;

while err >= tol || it >= 1e7

    r    = (b - A*x)
    x    = x + r

    err  = norm(x-xold);
    xold = x;
    it   = it + 1;

end

x 
% Running iteration by iteration shows that our guessing is too big!