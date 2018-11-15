function [ Q, R, L ] = Algorithm3( A )

Implemented in NREL
QR decomposition using two synch classical GS .
Sample use: A = rand(100); Algorithm3(A);

n=size(A, 1);
m=size(A, 2);


Q = zeros(n,m);
T = zeros(m,m);
R = eye(m,m);


for j=1: m
a = A(:, j);
Q(:, j) = a;

if j>1
    T(1:j-1, j-1) = Q(:, 1:j-1)'*Q(:, j-1);
    L(1:j-1, 1:j-1) = tril(T(1:j-1, 1:j-1),-1);
    R(1:j-1, j) = (eye(j-1,j-1) - L(:, 1:j-1)-L(:, 1:j-1)')*Q(:, 1:j-1)'*a;
    Q(:,j) = Q(:,j) - Q(:, 1:j-1)*R(1:j-1,j)-Q(:,1:j-1)*Q(:,1:j-1)'*(Q(:,j)-Q(:,1:j-1)*R(1:j-1,j ));
    
end
R( j, j) = norm(Q(:,j));
Q(:, j) =  Q(:, j)/norm(Q(:,j));




end

fprintf('||Q(:,1:%d)^TQ(:, 1:%d) - I|| = %16.16e|| \n', m, m, norm(Q(:, 1:m )'*Q(:, 1:m)-eye(m,m)));
 fprintf('||A - QR||/||A|| = %16.16e \n\n', norm(A-Q*R, 'fro')/norm(A, 'fro'));

end