% 101x +  12y - 13z  = 14
%  21x + 201y + 23z  = 24
% -31x +  32y + 301z = 34
%  Ax = b 
%  x = A^(-1)*b

% Assemble matrices and vectors
% Matrix A must be diagonally dominant! 
A = [101 12  -13  ;...
      21 201 23   ;...
     -31 32  301 ];

b = [ 14 ;
      24 ;
      34 ];

% Solve
% x = A^(-1)*b;
% or
% x = inv(A)*b;
% or
% x = linsolve(A,b);
% or
x = A\b;

% Check results
x = x' % transpose to be a row vector
Result1 = int8(sum(x .* A(1,:))) == b(1)
Result2 = int8(sum(x .* A(2,:))) == b(2)
Result3 = int8(sum(x .* A(3,:))) == b(3)
