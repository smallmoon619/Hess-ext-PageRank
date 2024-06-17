function [x,r,mv] = RestartArnoldi(A,d,n,alpha,k,v1,tol,p,mv)
%% run thick restarted arnoldi steps 2-7
%% run arnoldi process for H,V,H_bar
%r = 1;
iter = 0;
p0 = p;
%v1 = v;
while iter < 2
    if iter == 0
        [V,H,mv] = ArnoldiProcess(A,v1,k,n,d,alpha,mv);
    else
        [V,H,mv] = STAArnoldiProcess(A,v1,V,H,p,k,n,d,alpha,mv);
    end
    H_bar = H(1:k,:);
    p = p0;
    [eigvector,eigvalue] = eig(H_bar);
   [eigvalue,eigvector]=PTselectevr(eigvalue,eigvector,k,p);
   y1 = eigvector(:,1);
   x=V*(H*eigvector(:,1));
    oldx=V(:,1:k)*eigvector(:,1);
   r=norm(x-oldx,2);
    if H(k+1,k)*abs(y1(k))<tol
        x=V(:,1:k)*y1;
        x = x/sum(x);
        break;
    end
    clear y1;
%    oldx = x;
%    x=V(:,1:k)*eigvector(:,1);
   x = x/sum(x);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% varying p if necessary
   if imag(eigvalue(p))>0
       p=p-1;%decrease p, alternatively, we may increase p by 1 
       eigvalue=eigvalue(1:p);
       eigvector=eigvector(:,1:p);
   end
   iter=iter+1;
   if r<=tol
        x=V(:,1:k)*eigvector(:,1);%approximate eigenvectors corresponding to the largest eigenvalue
        x = x/sum(x);
        break
    end
    %% 划分向量的虚部与实部
    j=1;
    Wp=[];
    while j<=p
     if imag(eigvalue(j)) == 0%real
         Wp=[Wp eigvector(:,j)];
     end
     if imag(eigvalue(j)) ~= 0%complex
         Wp=[Wp real(eigvector(:,j))];
         j1=j+1;            
       if j1<=p
         if (real(eigvalue(j))==real(eigvalue(j1)))&(imag(eigvalue(j))==-imag(eigvalue(j1)))&imag(eigvalue(j))~=0
            Wp=[Wp imag(eigvector(:,j))];
            j=j1;               
         end
       end
     end   
     j=j+1;         
    end
    %% 正交化向量
    Wp = Gram_Schmidt_Orthogonalization(Wp);
    Wp1 = [Wp;zeros(1,p)];
    s=[zeros(k,1);1];
    Wp1 = [Wp1 s];
    %% 生成V^new,H^new
    V = V*Wp1;
    H = Wp1'*H*Wp;
    %% 返回step3 of the thick restart algorithm vp+1 = v
    v1 = V(:,p+1);
end