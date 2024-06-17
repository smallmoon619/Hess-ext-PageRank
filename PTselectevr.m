function [ev,evr]=PTselectevr(ev,evr,m,l)
%the function is to sort and choose l wanted Ritz pairs,i.e.,ev and evr
ev=diag(ev);
for i=2:m
 j=m-i+1;
   for k=1:j
       p=k+1;
       if abs(ev(k))<abs(ev(p))%da
         key=ev(k);
         ev(k)=ev(p);
         ev(p)=key;
         key=evr(:,k);
         evr(:,k)=evr(:,p);
         evr(:,p)=key;
      end
   end
end
for k=1:m-1
    p=k+1;
    x1=real(ev(k));
    x2=real(ev(p));
    y1=imag(ev(k));
    y2=imag(ev(p));
    if (x1==x2)&(y1==-y2)&(y1<0)
        ev(k)=conj(ev(k));
        ev(p)=conj(ev(p));
        evr(:,k)=conj(evr(:,k));
        evr(:,p)=conj(evr(:,p));
    end
end
ev=ev(1:l,1);
%ev=ev+shift*ones(l,1);%theta,将近似特征值还原
evr=evr(:,1:l);
%s=evr(:,1);

clear key;
clear x1;
clear x2;
clear y1;
clear y2;
clear p;