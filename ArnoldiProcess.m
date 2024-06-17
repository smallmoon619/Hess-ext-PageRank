function [V,H,mv]=ArnoldiProcess(A,v1,k,n,d,alpha,mv)
v = ones(n,1)/n;
V=zeros(n,k+1);
H=zeros(k+1,k);
e = ones(n,1);
V(:,1) = v1/norm(v1,2);
for j = 1:k
    z = alpha * (A' * V(:,j) + (d' * V(:,j)) * v) + (1 - alpha)*(e'*V(:,j))*v;
    mv = mv + 1;
    for i=1:j
        H(i,j) = z' * V(:,i);
        z = z - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(z,2);
    if H(j+1,j) ==0
        disp('break down');
        break;
    end
    if j <= k
        V(:,j+1) = z/H(j+1,j);
    end
end