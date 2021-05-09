function [res]= accuFucV2(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,numclass);
maxIter = 50;
tmp1 = zeros(maxIter,1);
tmp2 = zeros(maxIter,1);
for iter = 1:maxIter
    indx = litekmeans(U_normalized,numclass,'MaxIter',100, 'Replicates',10);
    indx = indx(:);
    [newIndx] = bestMap(Y,indx);
    tmp1(iter) = mean(Y==newIndx);
    tmp2(iter) = MutualInfo(Y,newIndx);
end
res = [max(tmp1),max(tmp2)];