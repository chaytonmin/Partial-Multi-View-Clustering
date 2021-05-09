function [H_normalized,theta,objective]= mkkmeans_train(Km, cluster_count,qnorm)

numker = size(Km, 3);
theta = ones(numker,1)/numker;
K_theta = mycombFun(Km, theta.^qnorm);

opt.disp = 0;
iteration_count = 0;
flag =1;
while flag
    iteration_count = iteration_count+1;
    fprintf(1, 'running iteration %d...\n', iteration_count);
    [H, ~] = eigs(K_theta, cluster_count, 'LA', opt);
%     Q = zeros(numker);
%     for m = 1:numker
%         Q(m, m) = trace(Km(:, :, m)) - trace(H' * Km(:, :, m) * H);
%     end
%     res = mskqpopt(Q, zeros(numker, 1), ones(1, numker), 1, 1, zeros(numker, 1), ones(numker, 1), [], 'minimize echo(0)');
%     theta = res.sol.itr.xx;
    [theta] = updateabsentkernelweightsV2(H,Km,qnorm);
    K_theta = mycombFun(Km, theta.^qnorm);
    objective(iteration_count) = -trace(H' * K_theta * H) + trace(K_theta);
    if iteration_count>2 && (abs((objective(iteration_count-1)-objective(iteration_count))/(objective(iteration_count-1)))<1e-3...
            || iteration_count>100)
        flag =0;
    end
end
% H_normalized = H ./ repmat(sqrt(sum(H.^2, 2)), 1, cluster_count);
H_normalized = H;