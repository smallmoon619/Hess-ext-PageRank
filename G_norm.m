function Gnorm = G_norm(r,u)
% r is a diagonal vector of diagonal matrix G, column vector
% u is a vector, column
% ��������G����
qx_sq = u.^2;
u_sum = r' * qx_sq;
Gnorm = sqrt(u_sum);