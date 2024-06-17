function [iter,mv,time,res] = Arnoldi_Inout(A,x,tol,eta,k,p,maxit,alpha,d,n)
%% input
% k the maximun size of the subspace 
% x inital vector
% tol outer tolerance
% eta inner tolerance
tic;
v = ones(n,1)/n; % positive vector                
alpha1 = alpha - 0.1; 
alpha2 = alpha - 0.1;
beta = 0.5;
r = 1;rin = 1; rin0 = rin; 
res = [];
iter = 0;
mv = 0;
%% Power-Inner-Outer_Arnoldi method
while r > tol
    restart=0;
    v1 = x;
    [x,residual,mv] = RestartArnoldi(A,d,n,alpha,k,v1,tol,p,mv);  
    if residual < tol
        break
    end
    while restart < maxit && r > tol
        x= x/norm(x,1);
        z = A' * x + (d' * x) * v;  mv = mv + 1;
        r = norm( alpha * z + (1 - alpha) * v-x,2)/norm(x,2);
        r0 = r;r1 = r;ratio = 0;
        while ratio < alpha1 && r > tol
            f = (alpha-beta)*z + (1-alpha)*v;
            ratio1=0;
            while ratio1 < alpha2 && rin > eta
                x = f+beta*z;
                z = A' * x + (d' * x) * v;  mv = mv + 1;
                rin = norm(f+beta*z-x,2)/norm(x,2);
                ratio1 = rin/rin0;
                rin0 = rin;
            end
            r = norm(alpha * z + (1 - alpha) * v-x,2)/norm(x,2);
            ratio = r/r0;
            r0 = r;
        end
         x = alpha * z + (1 - alpha) * v;
         res = [res;r];
         iter = iter + 1;
         x = x / norm(x,1);
         if r/r1 > alpha1
             restart = restart + 1;
         end
    end  
end
time = toc;