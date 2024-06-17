function [x,mv]= EXTprocedure (A,x0,lambda2,lambda3,mv,alpha,d,v) 
z = A' * x0 + (d' * x0) * v;
x1 = alpha * z + (1 - alpha) * v;
r1=norm(x1 - x0,2);
mv = mv + 1;
if imag(lambda2)==0 && imag(lambda3)==0
    z = A' * x1 + (d' * x1) * v;
    x2 = alpha * z + (1 - alpha) * v;
    mv = mv + 1;
    x = x2 -(lambda2+lambda3) * x1  + lambda2 * lambda3 * x0;
    r2 = norm(x-x2,2);
    if r1<r2
        x=x2; 
    end
elseif imag(lambda2+lambda3)==0
    z = A' * x1 + (d' * x1) * v;
    x2 = alpha * z + (1 - alpha) * v;
    mv = mv + 1;
    x = x2 -2 * real(lambda2) * x1  + abs(lambda2)^2 * x0;
    r2 = norm(x - x2,2);
    if r1<r2
        x=x2; 
    end
elseif imag(lambda2)==0 && imag(lambda3)~=0
    x = x1 - lambda2 *  x0;
    r2 = norm(x - x1,2);
    if r1<r2
        x=x1; 
    end
elseif imag(lambda2)~=0
    x = x1 - alpha *  x0;
    r2 = norm(x - x1,2);
    if r1<r2
        x=x1;
    end
end
x = x/norm(x,1);
end    