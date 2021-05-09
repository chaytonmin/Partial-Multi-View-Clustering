function [gamma]= updateabsentkernelweightsV2(T,K,qnorm)

num = size(K,1);
nbkernel = size(K,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0 = eye(num)-T*T';
a = zeros(nbkernel,1);
for p = 1 : nbkernel
    a(p) = trace( K(:,:,p) * U0);
end
gamma = a.^(-1/(qnorm-1))/sum(a.^(-1/(qnorm-1)));
gamma(gamma<eps)=0;
gamma = gamma/sum(gamma);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q = zeros(nbkernel);
% for p = 1:nbkernel
%     Q(p, p) = trace(K(:, :, p)) - trace(T' * K(:,:,p) * T);
% end
% res = mskqpopt(Q, zeros(nbkernel, 1), ones(1, nbkernel), 1, 1, zeros(nbkernel, 1), ones(nbkernel, 1), [], 'minimize echo(0)');
% gamma = res.sol.itr.xx;