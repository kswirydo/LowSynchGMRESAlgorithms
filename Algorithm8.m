
% this code is not restarted - not meant for large matrices (!!!)
% Example use: A = rand(1000,1); b = rand(1000,1); x=zeros(1000,1);
% Algorithm8(A,b,x0,1e-10);
% (c) Julien Langou (CU Denver), K. Swirydowicz (NREL), S. J. Thomas (NREL)

function [r,normhistory ] =  Agorithm8(A, b, x0, tol)
r = b-A*x0;
tol = tol * norm(r);
[n m] = size(A);
V = zeros(n,n);
H = zeros(n+1,n);

normhistory = zeros(n+2,1);
cs(1:n)      = zeros(n,1);
sn(1:n)      = zeros(n,1);
conv = 0;

V(:,1) = r/norm(r);
i=0;
e1 = zeros(n,1);
e1(1) = 1.0;
s      = norm(r)*e1;

normhistory(1) = norm(r);
while ((norm(r)>tol)&&conv == 0);
    i = i+1;
    z = A*V(:,i);
    w = z;
    
    for j=1:i
        H(j,i) = w'*V(:,j);
        w = w-H(j,i)*V(:, j);
    end
    
    H(i+1,i) = norm( w );
    if H(i+1,i) > 0
        V(:,i+1) = w / H(i+1,i);
        
        % Apply Givens rotation
        for k = 1:i-1
            temp     =  cs(k)*H(k,i) + conj(sn(k))*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + conj(cs(k))*H(k+1,i);
            H(k,i)   = temp;
        end
        
        % Form i-th rotation matrix
        [ cs(i), sn(i) ] = rotmat( H(i,i), H(i+1,i) );
        
        % Approximate residual norm
        temp        = cs(i)*s(i);
        s(i+1)      = -sn(i)*s(i);
        s(i)        = temp;
        H(i,i)      = cs(i)*H(i,i) + conj(sn(i))*H(i+1,i);
        H(i+1,i)    = 0.0;
        normhistory(i+1) = abs(s(i+1));
        
        
        
        if normhistory(i) <= tol
            conv = 1;
            
        end
    else
        conv = 1;
    end
    
end
i=i-1;
y = H(1:i,1:i) \ s(1:i);
x = x0 + V(:,1:i)*y;
r = b-A*x;


normhistory(i+2) = norm(r);
normhistory = normhistory(1:i+2);

close all
figure(1)
plot(log10(normhistory), 'b*-');
xlim([1 i+2]);
xlabel('iteration');
ylabel('log_{10} of the residual');
grid on;

fprintf("Converged in %d iters, norm of the final residual %1.16e \n", i, norm(r));


end
