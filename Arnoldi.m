function [iter,mv,time,res] = Arnoldi(A,x,tol,k,alpha,d,n)
% x initial vector
% Arnoldi process with the 2-norm
%%
tic
iter = 0;
mv = 0;
I = [eye(k,k);zeros(1,k)];
%-----------------------------------
% v = ones(n,1)/n;
% prev_x = x;
% z = A' * x + (d' * x) * v;
% x = alpha * z + (1 - alpha) * v;
% residual = norm(x-prev_x,1)/norm(x,1);
% mv = mv + 1;iter = iter + 1;
% res = [residual];
res = [];
residual = 1;
%-----------------------------------
while residual > tol
    [Q,H,mv]=ArnoldiProcess(A,x,k,n,d,alpha,mv);% HÊÇ(k+1) x k,QÊÇn x (k+1)
    [U,S,V]=svd(H-I);
    x = Q(:,1:k) * V(:,k);
    r = S(k,k) * Q * U(:,k);
%     residual = norm(r,1)/norm(x,1);
%     residual = norm(r,2);
    residual = norm(r,2)/norm(x,2);
    iter = iter + 1;
    res = [res; residual];
end
time = toc;

% residual = S(k,k);
%     r = S(k,k) * Q * U(:,k);
%     residual = norm(r,1);
    %residual = norm(x-prev_x,1);
