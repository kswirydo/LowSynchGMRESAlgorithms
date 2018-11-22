function [ Q, R ] = Algorithm1( A, iter )

% Implemented in NREL
% QR decomposition using Ruhe notation classical Gram-Schmidt
% Sample use: A = rand(100); Algorithm1(A. 2);
% second parameter is the number of orth steps per vector
% (c) Julien Langou (CU Denver), K. Swirydowicz (NREL), S. J. Thomas (NREL)

n=size(A, 1);
m=size(A, 2);


Q = zeros(n,m);
R = eye(m,m);
R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/norm(A(:,1));

for i=2: m
    a = A(:, i);
    R(1:i-1, i) = zeros(i-1,1);
    for k=1:iter
        s = Q(:, 1:i-1)'*a;
        R(1:i-1, i) = R(1:i-1, i) +s;
        a = a - Q(:, 1:i-1)*s;
        
    end
   
    R( i, i) = norm(a);
    Q(:, i) =  a/R(i,i);
   
end

fprintf('||Q(:,1:%d)^TQ(:, 1:%d) - I|| = %16.16e|| \n', m, m, norm(Q(:, 1:m )'*Q(:, 1:m)-eye(m,m)));
fprintf('||A - QR||/||A|| = %16.16e \n\n', norm(A-Q*R, 'fro')/norm(A, 'fro'));

end
