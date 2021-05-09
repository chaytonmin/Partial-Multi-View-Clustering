function [Kyr]= absentKernelImputation(Kx,Kycc,mset,alpha0)

%% missing the indices of missing samples.
n = size(Kx,1);
n0 = length(mset);
%% c is the indices of available samples.
cset = setdiff(1:n,mset);
%%
Kx0 = Kx([cset,mset],[cset,mset]);
clear Ky
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dx = diag(sum(Kx0));
% Lx = eye(n) - Dx\Kx0;
Lxmm = Kx0(n-n0+1:end,n-n0+1:end);
Lxmm = (Lxmm+Lxmm')/2;
Lxcm = Kx0(1:n-n0,n-n0+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lxcmmm = -Lxcm/(Lxmm + alpha0*eye(n0));
Kycm = Kycc*Lxcmmm;
Kymm = Lxcmmm'*Kycm;
Kyr0 = [Kycc Kycm; Kycm' Kymm];
Kyr0 = (Kyr0+Kyr0')/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[val,indxxx] = sort([cset,mset],'ascend');
Kyr=Kyr0(indxxx,indxxx)+1e-12*eye(n);
% Kyr = Kyr/trace(Kyr);