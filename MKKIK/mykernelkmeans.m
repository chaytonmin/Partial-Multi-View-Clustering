function [H_normalized,obj]= mykernelkmeans(K,cluster_count)

K = (K+K')/2;
opt.disp = 0;
[H,~] = eigs(K,cluster_count,'la',opt);
obj = trace(H' * K * H) - trace(K);
% H_normalized = H ./ repmat(sqrt(sum(H.^2, 2)), 1,cluster_count);
H_normalized = H;