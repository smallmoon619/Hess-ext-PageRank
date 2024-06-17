clear all;
% clc;
%% %%%导入数据，数据是mat格式%%%%%%
% load('soc-Slashdot0902.mat')
load('ljournal-2008.mat')
% load('amazon0505.mat')    
% load('wiki-Talk.mat')   
% load('wikipedia-20070206.mat') 
A = Problem.A;
clear Problem;
%% A为超链接矩阵即transition matrix%%%%%
A = normalwebmatrix(A);
%% %%%%%%%%%%寻找悬挂节点,若节点i为悬挂节点，则d(i)=1,否则为0%%%%%%%%
n = max(size(A));
rowsumvector = ones(1,n)*A';%行和向量
nonzerorows = find(rowsumvector);%寻找行和不为零的行下标
zerorows = setdiff(1:n,nonzerorows); %返回从1到n里nonzerorow没有的数
num_d = length(zerorows);%悬挂节点的个数
d = sparse(zerorows,ones(num_d,1),ones(num_d,1),n,1);%悬挂向量
%%
eta = 1e-2;
tol = 1e-8;
x = ones(n,1)/n;
k = 10;% the dimension of Krylov subspace
% l=40;
maxit = 8;% The number of extrapolation procedure
maxit1 = 4;% Control parameter of Arnoldi_Inout
p = 2;% Control parameter of Arnoldi_Inout
%%
alpha = 0.85;
fprintf('%f \n',alpha)
[iter_1,mv_1,time_1,res_1] = Arnoldi(A,x,tol,k,alpha,d,n); % the Arnoldi-type algorithm
[iter_2,mv_2,time_2,res_2] = Arnoldi_Inout(A,x,tol,eta,k,p,maxit1,alpha,d,n);% the Arnoldi-Inout algorithm
[iter_5,mv_5,time_5,res_5] = A_Arnoldi(A,x,tol,k,alpha,d,n);%  the adaptively accelerated Arnoldi algorithm
[iter_6,mv_6,time_6,res_6] = Power_GArnoldi(A,x,tol,k,maxit1,alpha,d,n);%the adaptive Power-GArnoldi algorithm
[iter_7,mv_7,time_7,res_7] = HessenPR(A,x,k,n,d,alpha,tol);% the Hessenberg-type algorithm
[iter_8,mv_8,time_8,res_8] = Hessen_Ex1(A,x,tol,k,alpha,d,n,maxit); % the Hessenberg-extrapolation algorithm
fprintf('%d(%d) %d(%d) %d(%d) %d(%d) %d(%d) %d(%d) \n',iter_1,mv_1,iter_2,mv_2,iter_5,mv_5,iter_6,mv_6,iter_7,mv_7,iter_8,mv_8)
fprintf('%f %f %f %f %f %f \n',time_1,time_2,time_5,time_6,time_7,time_8)
res_1 = [1;res_1];res_2 = [1;res_2];res_5 = [1;res_5];res_6 = [1;res_6];res_7 = [1;res_7];res_8 = [1;res_8];
subplot(2,2,1)
semilogy([0:1:length(res_1)-1],res_1,'k-.','linewidth',1.8)
hold on
semilogy([0:1:length(res_2)-1],res_2,'m-.','linewidth',1.8)
semilogy([0:1:length(res_5)-1],res_5,'g-.','linewidth',1.8)
semilogy([0:1:length(res_6)-1],res_6,'b-.','linewidth',1.8)
semilogy([0:1:length(res_7)-1],res_7,'y-.','linewidth',1.8)
semilogy([0:1:length(res_8)-1],res_8,'r-.','linewidth',1.8)
ylim([1e-10 1e0])
title('\alpha=0.85')
xlabel('Iterations')
xlabel('Iterations')
ylabel('Residual norms')
legend('Arnoldi','AIO','A-A','A-PGA','Hessenberg','Hess-Ext')

% print amazon0505.eps -depsc2 -r800
% print(gcf,'-dpng','abc.png') 