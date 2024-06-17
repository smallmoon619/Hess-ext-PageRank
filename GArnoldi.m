function [x,residual,mv] = GArnoldi(A,r,x,k,n,d,alpha,mv)
%% 利用
for i = 1:2
    I = [eye(k,k);zeros(1,k)];
    [Q,H,mv] = G_ArnoldiProcess(A,r,x,k,n,d,alpha,mv);  % H是(k+1) x k,Q是n x (k+1)
    [U,S,V]=svd(H-I);
    x = Q(:,1:k) * V(:,k);
    r = S(k,k) * Q * U(:,k);
    residual = norm(r,2);
    r = abs(r) / norm(r,1);
end