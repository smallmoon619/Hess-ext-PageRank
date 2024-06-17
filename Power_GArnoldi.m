function [iter,mv,time,res] = Power_GArnoldi(A,x,tol,k,maxit,alpha,d,n)
% x initial vector
% r residual vector
% Arnoldi process with the G-norm
v = ones(n,1)/n;
tau = ones(n,1);
% x_G_norm = G_norm(tau,x);
% x = x/x_G_norm;
beta = alpha - 0.1;
r = 1;
%e = ones(n,1);
%% run arnoldi process for H,V,H_bar
res = [];
iter = 0;
mv = 0;
%% Power_GArnoldi method
tic;
while r > tol
    restart=0; 
    [x,residual,mv] = GArnoldi(A,tau,x,k,n,d,alpha,mv);
    if residual < tol
        disp('convergence')
        x = x/norm(x,1);
        break
    end
    while restart < maxit && r > tol
        x=abs(x)/norm(x,1);
        r0 = r;
        ratio = 0;
        r1 = r;
        while ratio < beta && r > tol
            prev_x = x;
            z = A' * x + (d' * x) * v;
            mv = mv + 1;
            x = alpha * z + (1 - alpha) * v;
            tau = x-prev_x;
            r = norm(tau,2)/norm(prev_x,2);
            res=[res; r];
            iter = iter + 1;
            if r <= tol
                break;
            end
            ratio = r/r0;
            r0 = r;
        end
        if r/r1 > beta
            restart = restart + 1;
        end
    end  
    tau = abs(tau)/norm(tau,1);
end
% res = [1;res];
time = toc;