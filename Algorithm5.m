function [ Q, R, L ] = Algorithm5( A )

% Implemented in NREL
% QR decomposition with triangular solve lagged mod GS .
% Sample use: A = rand(100); Algorithm5(A);
% (c) Julien Langou (CU Denver), K. Swirydowicz (NREL), S. J. Thomas (NREL)

n=size(A, 1);
m=size(A, 2);


Q = zeros(n,m);
T = zeros(m,m);
R = eye(m,m);

for j=1: m
    a = A(:, j);
    Q(:, j) = a;
    if j>1
           
        tmp = Q(1:n,1:j-1)'*Q(1:n,j-1:j); %%%%%%%%% <- this is the only synchronization
        %
        T(1:j-2,j-1) = tmp(1:j-2,1);
        R(j-1,j-1) = tmp(j-1,1);
        R(1:j-1,j) = tmp(1:j-1,2);
        
        R(j-1,j-1) = sqrt( R(j-1,j-1) );
        Q(1:n,j-1) = Q(1:n, j-1)/R(j-1,j-1);
        
        R(j-1,j) = R(j-1,j) / R(j-1,j-1);
        if j>2, T(1:j-2, j-1) = T(1:j-2, j-1)/R(j-1,j-1); end
        
        T(j-1,j-1) = 1.0;
        R(1:j-1,j) = T(1:j-1, 1:j-1)'\R(1:j-1,j);
        
        Q(1:n,j) = Q(1:n,j) - Q(1:n,1:j-1)*R(1:j-1,j);
        
    end
   
end
R(m,m) = norm(Q(:,m));
Q(:,m) = Q(:,m)/R(m,m);

fprintf('||Q(:,1:%d)^TQ(:, 1:%d) - I|| = %16.16e|| \n', m, m, norm(Q(:, 1:m )'*Q(:, 1:m)-eye(m,m)));
fprintf('||A - QR||/||A|| = %16.16e \n\n', norm(A-Q*R, 'fro')/norm(A, 'fro'));

end
