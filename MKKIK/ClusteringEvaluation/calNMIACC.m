function [acc1,nmi1]= calNMIACC(K,Y,numclass)

maxIter = 10;
tmpnmi = zeros(maxIter,1);
tmpacc = zeros(maxIter,1);
dis = zeros(maxIter,1);
for iter = 1:maxIter
    [U] =  mykernelkmeans(K,numclass);
    [indx,C,SUMD] = kmeans(U,numclass);
    indx = indx(:);
    [newIndx] = bestMap(Y,indx);
    tmpnmi(iter) = MutualInfo(Y,newIndx);
    tmpacc(iter) = mean(Y==newIndx);
    dis(iter) = sum(SUMD);
end
nmi1 = max(tmpnmi);
acc1 = max(tmpacc);