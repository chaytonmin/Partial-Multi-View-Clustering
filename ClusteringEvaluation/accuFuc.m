function [acc]= accuFuc(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,numclass);

% maxIter = 50;
% tmp = zeros(maxIter,1);
% for iter = 1:maxIter
%     %% [indx] = knkmeans(K,numclass);
%     indx = litekmeans(U_normalized,numclass, 'MaxIter',100, 'Replicates',10);
%     indx = indx(:);
%     [newIndx] = bestMap(Y,indx);
%     tmp(iter) = mean(Y==newIndx);
% end
% % acc = [mean(tmp);std(tmp)];
% acc = max(tmp);
% 

indx = litekmeans(U_normalized,numclass, 'MaxIter', 100, 'Replicates',50);
indx = indx(:);
[newIndx] = bestMap(Y,indx);
acc = mean(Y==newIndx);