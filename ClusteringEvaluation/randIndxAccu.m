function [accu0]= randIndxAccu(X,Y)

maxIter = 10;
tmp1 = zeros(maxIter,1);
kcluster = length(unique(Y));
num1 = length(Y);
for iter = 1:maxIter
    opts = statset('Display','final');
    indx = kmeans(X,kcluster,'EmptyAction','singleton', 'Replicates',10, 'Options',opts);
    account0 = 0;
    for i = 1: num1
        for j = i+1:num1
            if (Y(i)==Y(j) && indx(i)==indx(j)) || (Y(i)~=Y(j) && indx(i)~=indx(j))
                account0 = account0 +1;
            end
        end
    end
    tmp1(iter) = 2*account0/(num1*(num1-1));
end
accu0 = max(tmp1);
