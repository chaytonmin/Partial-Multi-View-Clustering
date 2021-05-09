function [result] = MultiNMF_incomplete_original_l21(X, C, K, label, options)
% This is a module of Multi-View Non-negative Matrix Factorization(MultiNMF)
%
% Notation:
% X ... a cell array containing all views for the data. Each
% C.... 对角权重矩阵W，
% K ... number of hidden factors
% label ... ground truth labels
% options ... a cell containing the parameters
% Depend on the dataset and the normalization of the dataset, the two parameters alpha and beta need to be tuned. 
% implemented by Weixiang Shao (wshao4@uic.edu)

viewNum = length(X);
Rounds = options.rounds;

%%将W转化为需要的对角形式
for i=1:viewNum
    C{i}=diag(C{i});
end

U_ = [];
V_ = [];

U = cell(1, viewNum); 
V = cell(1, viewNum);

j = 0;
log = 0;
ac = 0;

while j < 3%原来是3
    j = j + 1;
    if j == 1
        [V{1}, U{1}] = NMF(X{1}', K, options, V_, U_);
    else
        [V{1}, U{1}] = NMF(X{1}', K, options, V_, U{viewNum});
    end
    printResult(U{1}, label, K, options.kmeans);
    for i = 2:viewNum
        [V{i}, U{i}] = NMF(X{i}', K, options, V_, U{i-1});
        printResult(U{i}, label, K, options.kmeans);
    end
end

optionsForPerViewNMF = options;
oldL = 10000000;
oldU = U;
oldV = V;
tic
j = 0;
oldcentroidU = zeros(size(C{1},2), K);
converge = 0;
while j < Rounds
    j = j + 1;
        %更新U*
        CU = options.alpha(1)*(C{1}.^2)*U{1};%求u*第一项
        CC = options.alpha(1)*(C{1}.^2);%求U*第二项

        for i = 2:viewNum
            CU = CU + options.alpha(i)*(C{i}.^2)*U{i};
            CC = CC + options.alpha(i)*(C{i}.^2);
        end
        CC_inv = diag(1./diag(CC));
        centroidU = CC_inv*CU;%U*
    logL = 0;
    for i = 1:viewNum
        tmp1 = C{i}*(X{i} - U{i}*V{i}');
        tmp2 = C{i}*(U{i} - centroidU);
        tmp3 = 0;
        for k =1:size(U{i},2);
            tmp3 = tmp3 + norm(U{i}(:,k));
        end
        logL = logL + sum(sum(tmp1.^2)) + options.alpha(i) * sum(sum(tmp2.^2)) + options.beta(i)*tmp3;
    end
    log(end+1) = logL;
    logL;
    if(oldL < logL)
        j = j;
        disp('objective function value increasing');
        ac(end+1) = printResult(centroidU, label, K, options.kmeans);%(end+1)表示下标自动扩大
    else
        ac(end+1) = printResult(centroidU, label, K, options.kmeans);
    end
    avg_diff = sum(sum(abs(oldcentroidU - centroidU).^2));
    fprintf('The average diff is %d for iteration %f\n', avg_diff, j);
    if  avg_diff< 1e-12
        converge = 1;
        fprintf('converge at iteration %f\n', j);
    end
    oldU = U;
    oldV = V;
    oldL = logL;
    oldcentroidU = centroidU;
    %更新U,V
    for i = 1:viewNum
        optionsForPerViewNMF.alpha = options.alpha(i);
        optionsForPerViewNMF.beta = options.beta(i);
        [U{i}, V{i}] = PerViewNMF_incomplete_original_l21(X{i}, K, centroidU, optionsForPerViewNMF, U{i}, V{i}, C{i});
    end
end
%计算最终的评价结果
[result] = ComputeStats(centroidU, label, K,options.kmeans);
toc
