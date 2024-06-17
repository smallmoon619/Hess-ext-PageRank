function [V,H,mv] = G_ArnoldiProcess(A,r,x,k,n,d,alpha,mv)
% the initial matrix G = I = diag(r)
% the code for the paper of YiGuoJian-2012, Method 3.1
% Arnoldi process with the G-norm
v = ones(n,1)/n;
V=zeros(n,k+1);
H=zeros(k+1,k);
e = ones(n,1);  
x_G_norm = G_norm(r,x);
V(:,1) = x/x_G_norm;
for j = 1:k
    z = alpha * (A' * V(:,j) + (d' * V(:,j)) * v) + (1 - alpha)*(e'*V(:,j))*v;
    mv = mv + 1;
    for i=1:j
        qk = r.*V(:,i);
        H(i,j) = z' * qk;
        z = z - H(i,j)*V(:,i);
    end
    H(j+1,j) = G_norm(r,z);
    if H(j+1,j) ==0
        disp('break down');
        break;
    end
    if j <= k
        V(:,j+1) = z/H(j+1,j);
    end
end
