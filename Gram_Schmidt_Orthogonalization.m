function T = Gram_Schmidt_Orthogonalization(p_tr)
% һ��Ϊһ������
[row,col]= size(p_tr);
T = zeros(row,col);
T(:,1)=p_tr(:,1);
for i = 2 : col
    for j = 1: i-1
        p_tr(:,i)= p_tr(:,i) - ((T(:,j)' * p_tr(:,i))/(T(:,j)' * T(:,j))) * T(:,j);
    end
    T(:,i)=p_tr(:,i);
end
% ������λ��
for i = 1: col
    length=norm(T(:,i));
    for j = 1: row
        T(j,i)= T(j,i)/ length;
    end
end
end
