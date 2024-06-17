function [iter,mv,time,res] = Hessen_Ex1(A,x,tol,k,alpha,d,n,maxit)
tic;
mv = 0;
v = ones(n,1)/n;
e = ones(n,1);
iter = 0;
% p = zeros(n,1);
% p = 1:n;
%-----------------
L = zeros(n,k+1);
H = zeros(k + 1,k);
%-----------------
It = zeros(k + 1,k);
It(1:k,1:k) = eye(k);
%-----------------
r = 1;
res = [];
%-----------------
while r > tol
    p = 1:n;
    [~,i0] = max(abs(x)); 
    beta = x(i0); 
    L(:,1) = x/beta; 
    p(1) = i0;
    tmp = p(1);
    p(1) = p(i0);
    p(i0) = tmp;
    for j = 1:k
         u = alpha * (A' * L(:,j) + (d' * L(:,j)) * v) + (1 - alpha)*(e'*L(:,j))*v;
         mv = mv + 1;
        for i = 1 : j
             H(i,j) = u(p(i));
             u = u - H(i,j)*L(:,i);
        end 
%         if (j < n && ~isequal(u,zeros(n,1)))
%             [~,i0] = max(abs(u)); H(j+1,j) = u(i0);L(:,j + 1) = u/H(j + 1,j); p(j+1) = i0;
%         else
%             H(j+1,j) = 0;
%             break
%         end
        tp = zeros(n-j,1);
        if (j < n && ~isequal(u,zeros(n,1)))
            for i = j + 1:n
                tp(i-j) = u(p(i)); % Vector of u(w(k+1)),...,u(w(n))
            end
            [~,ind_p]= max(abs(tp)); % Maximum value of the vector abs(tp)
            i0 = j + ind_p;
            H(j + 1,j) = u(p(i0));
            L(:,j + 1) = u/H(j + 1,j);
            % Swap contents
            tmp = p(j+1);
            p(j + 1) = p(i0);
            p(i0) = tmp;
        else
            H(j+1,j) = 0;
        end
    end
    H1 = H - It;   % H_{m+1,m} - [I_m;0] = U*\Sigma*(V^T);
    [U,S,V] = svd(H1);
    x = (L(:,1:k))*V(:,k);
    deta = diag(S);
    deta = deta(end);
    qt = deta*(L*U(:,k)); 
    iter = iter + 1;
    r = norm(qt,2)/norm(x,2);
    res = [res;r];
    if r < tol
       break;
    end
    x = abs(x)/norm(x,1);
    H_bar = H(1:k,:);
    hm =sort( eig(H_bar),'descend');
    lambda2 = hm(2);
    lambda3 = hm(3);
    for i=1:maxit
       [x,mv]= EXTprocedure (A,x,lambda2,lambda3,mv,alpha,d,v);
    end 
    x = x/norm(x,1);
end
time = toc;
end