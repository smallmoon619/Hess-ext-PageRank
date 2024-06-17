function A = normalwebmatrix(G)
%% c = outdegree,r = indegree
[~,n] = size(G);
c = sum(G,1);    % each colunm sum
r = sum(G,2);    % each row sum
%% scale row/column to be 1 or 0,0 refers to no out links in this node
k = find(r ~= 0);
D = sparse(k,k,1./r(k),n,n);
%% obtain a matrix whose row/column sum is one or zero
% row stochastic matrix with some zero_row_sum
A = D * G;