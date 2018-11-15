function [ Q, R ] = Algorithm2( A, iter )

% Implemented in NREL
% QR decomposition using Ruhe notation mod Gram-Schmidt
% Sample use: A = rand(100); Algorithm2(A);


n=size(A, 1);
m=size(A, 2);


Q = zeros(n,m);
R = eye(m,m);
R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/norm(A(:,1));

for j=2: m
    a = A(:, j);
    R(1:j-1, j) = zeros(j-1,1);
    for i=1:j-1
        s = Q(:, i)'*a;
        a = a - s*Q(:,i);
        R(i, j) = R(i, j) +s;        
    end
  
    R( j, j) = norm(a);
    Q(:, j) =  a/R(j,j);

end

fprintf('||Q(:,1:%d)^TQ(:, 1:%d) - I|| = %16.16e|| \n', m, m, norm(Q(:, 1:m )'*Q(:, 1:m)-eye(m,m)));
fprintf('||A - QR||/||A|| = %16.16e \n\n', norm(A-Q*R, 'fro')/norm(A, 'fro'));

end