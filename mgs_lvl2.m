function [Q1, R1, T1] = mgs_lvl2(Q,R, T, j)
%this is the inner part of the loop in Alg 6
n = size(Q, 1);
if j>1
    
    
    tmp = Q(1:n,1:j-1)'*Q(1:n,j-1:j); %%%%%%%%% <- this is the only synchronization
    
    
    
    T(1:j-2,j-1) = tmp(1:j-2,1);
    R(j-1,j-1) = tmp(j-1,1);
    R(1:j-1,j) = tmp(1:j-1,2);
    
    R(j-1,j-1) = sqrt( R(j-1,j-1) );
    Q(1:n,j-1) = Q(1:n, j-1)/R(j-1,j-1);
    
    R(j-1,j) = R(j-1,j) / R(j-1,j-1);
    if j>2, T(1:j-2, j-1) = T(1:j-2, j-1)/R(j-1,j-1); end
    
    T(j-1,j-1) = 1.0;
    T(1:j-2, j-1) = (-1.0)*T(1:j-2, 1:j-2)*T(1:j-2, j-1);
    R(1:j-1,j) = T(1:j-1, 1:j-1)'*R(1:j-1,j);
    
    Q(1:n,j) = Q(1:n,j) - Q(1:n,1:j-1)*R(1:j-1,j);
    
end
Q1 = Q(:,1: j);
R1 = R(1:j,1:j);
T1 = T(1:j, 1:j);
end