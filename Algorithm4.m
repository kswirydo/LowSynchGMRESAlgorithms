function [ Q, R, L ] = Algorithm4( A )

% Implemented in NREL
% QR decomposition using mGS with trianglular solve insted of iteration.
% Sample use: A = rand(100); Algorithm4(A);

n=size(A, 1);
m=size(A, 2);


Q = zeros(n,m);
L = zeros(m,m);
R = eye(m,m);


for i=1: m
a = A(:, i);
Q(:, i) = a;

if i>1
    b = Q(:, 1:i-1)'*Q(:, i-1);
    L(1:i-1, i-1) = b; 
    
    R( 1:i-1, i) = (L(1:i-1, 1:i-1)')\(Q(:,1:i-1)'*a) ;
    Q(:, i) = Q(:, i) - Q(:, 1:i-1)*R( 1:i-1, i) ;
end
R( i, i) = norm(Q(:,i));
Q(:, i) =  Q(:, i)/norm(Q(:,i));




end

fprintf('||Q(:,1:%d)^TQ(:, 1:%d) - I|| = %16.16e|| \n', m, m, norm(Q(:, 1:m )'*Q(:, 1:m)-eye(m,m)));
 fprintf('||A - QR||/||A|| = %16.16e \n\n', norm(A-Q*R, 'fro')/norm(A, 'fro'));

end