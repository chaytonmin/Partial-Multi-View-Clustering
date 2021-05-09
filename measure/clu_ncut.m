function [idx,DistC,U,S] = clu_ncut(W,K)
% this routine groups the data X into K subspaces by NCut
% inputs:
%       W -- an N*N affinity matrix, N is the number of data points
%       K -- the number of subpaces (i.e., clusters)

n = size(W,1);
W(1:n+1:end) = 0;
W = (W + W')/2;
D = diag(1./(sqrt(sum(W,2))+eps));
LMat = eye(size(W))-D*W*D;
opts.tol = 1e-8;
opts.maxit = 30000;
[V,S]=eigs(LMat,K,'sr',opts);

U = V(:,1:K);
U = D*U;

[idx,~,~,DistC] = kmeans(U,K,'emptyaction','singleton','replicates',100,'display','off');
idx = idx';
