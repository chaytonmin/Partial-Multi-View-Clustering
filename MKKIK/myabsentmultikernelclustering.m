function [H_normalized,gamma,obj,KA] = myabsentmultikernelclustering(K,S,cluster_count,qnorm,algorithm_choose)

num = size(K,1);
nbkernel = size(K,3);
alpha0 = 1e-3;
%% S: num0*m, each column indicates the indices of absent samples
%% initialize kernel weights
gamma = ones(nbkernel,1)/nbkernel; 
%% initialize base kernels with zeros
if strcmp(algorithm_choose,'algorithm0')
    KA = feval(algorithm_choose,K,S);
else
    KA = feval(algorithm_choose,K,S);
end
%% combining the base kernels
KC  = mycombFun(KA,gamma.^qnorm);
flag = 1;
iter = 0;
while flag
    iter = iter + 1;
    %% update H with KC
%     fprintf(1, 'running iteration of the proposed algorithm %d...\n', iter);
    H = mykernelkmeans(KC,cluster_count);
   %% updata base kernels
    KA = zeros(num,num,nbkernel);
    for p =1:nbkernel
        if isempty(S{p}.index)
            KA(:,:,p) = K(:,:,p);
        else
            Kx = eye(num) - H*H';
            %             indexp = setdiff(1:nbkernel,p);
            %             for q = 1:nbkernel-1
            %                 Kx = Kx + lambda*gamma(indexp(q))*(-KA0(:,:,indexp(q)));
            %             end
            mis_indx = S{p}.index;
            obs_indx = setdiff(1:num,mis_indx);
            KA(:,:,p) = absentKernelImputation(Kx,K(obs_indx,obs_indx,p),mis_indx,alpha0);
        end
    end
    %% update kernel weights
    [gamma] = updateabsentkernelweightsV2(H,KA,qnorm);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [obj(iter)] = calObjV2(H,KA,gamma);
    %% KC  = mycombFun(KA,gamma.^qnorm);
    KC  = mycombFun(KA,gamma.^qnorm);
    if iter>2 && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-4 ||iter>30)
        flag =0;
    end
end
H_normalized = H./ repmat(sqrt(sum(H.^2, 2)), 1,cluster_count);