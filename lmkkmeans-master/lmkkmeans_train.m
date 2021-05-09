function [H_normalized,objective] = lmkkmeans_train(Km, cluster_count)

N = size(Km, 2);
P = size(Km, 3);
Theta = ones(N, P) / P;
K_Theta = calculate_localized_kernel_theta(Km, Theta);

opt.disp = 0;
iteration_count = 0;
flag =1;
while flag
    iteration_count = iteration_count +1;
    fprintf(1, 'running iteration %d...\n', iteration_count);
    [H, ~] = eigs(K_Theta, cluster_count, 'la', opt);
    HHT = H * H';
    
    Q = zeros(N * P, N * P);
    for m = 1:P
        start_index = (m - 1) * N + 1;
        end_index = m * N;
        Q(start_index:end_index, start_index:end_index) = eye(N, N) .* Km(:, :, m) - HHT .* Km(:, :, m);
    end
    res = mskqpopt(Q, zeros(N * P, 1), repmat(eye(N, N), 1, P), ones(N, 1), ones(N, 1), zeros(N * P, 1), ones(N * P, 1), [],...
        'minimize echo(0)');
    Theta = reshape(res.sol.itr.xx, N, P);
    K_Theta = calculate_localized_kernel_theta(Km, Theta);
    
    objective(iteration_count) = trace(H' * K_Theta * H) - trace(K_Theta);
    
    if iteration_count>2 && (abs((objective(iteration_count)-objective(iteration_count-1))/(objective(iteration_count)))<1e-3)
        flag =0;
    end
end
H_normalized = H ./ repmat(sqrt(sum(H.^2, 2)), 1, cluster_count);

%     stream = RandStream.getGlobalStream;
%     reset(stream);
%     state.clustering = kmeans(H_normalized, parameters.cluster_count, 'MaxIter', 1000, 'Replicates', 10);
%     state.objective = objective;
%     state.parameters = parameters;
%     state.Theta = Theta;
%     state.time = toc;