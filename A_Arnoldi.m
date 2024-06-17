function [it,mv,time,res] = A_Arnoldi(A,x,tol,k,alpha,d,n)
% x initial vector
% r residual vector
% the code for the paper of YiGuoJian-2012, Method 4.1 Adaptively Accelerated Arnoldi method 
% Arnoldi process with the G-norm
%%
r = ones(n,1);
it = 0;
mv = 0;
residual = 1;
res = [];
% x_G_norm = G_norm(r,x);
% x = x/x_G_norm;
tic
I = [eye(k,k);zeros(1,k)];
while residual > tol
    [Q,H,mv] = G_ArnoldiProcess(A,r,x,k,n,d,alpha,mv);% HÊÇ(k+1) x k,QÊÇn x (k+1)
    [U,S,V]=svd(H-I);
    x = Q(:,1:k) * V(:,k);
    r = S(k,k) * Q * U(:,k);
    residual = norm(r,2)/norm(x,2);
    it = it + 1;
    res = [res;residual];
    if residual < tol
        break;
    end
    r = abs(r)/norm(r,1);
end
time = toc;