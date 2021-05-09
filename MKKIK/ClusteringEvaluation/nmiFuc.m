function [acc]= nmiFuc(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
% maxIter = 20;
% U_normalized = U./ repmat(sqrt(sum(U.^2, 2)), 1,numclass);
% tmp = zeros(maxIter,1);
% for iter = 1:maxIter
%     %% kmeans(H_normalized, parameters.cluster_count, 'MaxIter', 1000, 'Replicates', 10);
%     indx = litekmeans(U_normalized,numclass,'MaxIter',100);
%     indx = indx(:);
%     [newIndx] = bestMap(Y,indx);
%     tmp(iter) = MutualInfo(Y,newIndx);
% end
% acc = [mean(tmp);std(tmp)];

[indx]= litekmeans(U,numclass, 'MaxIter',100, 'Replicates',30);
indx = indx(:);
[newIndx] = bestMap(Y,indx);
acc = MutualInfo(Y,newIndx);