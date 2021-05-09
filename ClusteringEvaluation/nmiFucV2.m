function [acc]= nmiFucV2(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
maxIter = 50;
U_normalized = U./ repmat(sqrt(sum(U.^2, 2)), 1,numclass);
tmp = zeros(maxIter,1);
for iter = 1:maxIter
    %% kmeans(H_normalized, parameters.cluster_count, 'MaxIter', 1000, 'Replicates', 10);
    indx = litekmeans(U_normalized,numclass,'MaxIter',100,'Replicates', 10);
    indx = indx(:);
    [newIndx] = bestMap(Y,indx);
    tmp(iter) = MutualInfo(Y,newIndx);
end
acc = max(tmp);