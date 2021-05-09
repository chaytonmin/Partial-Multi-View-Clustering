function [Ux Uy P2 P1 P3 F P R nmi avgent AR] = IMGclust(X2,Y2,X1,Y3,numClust,truth,option)

% Partial codes are from 
% S.-Y. Li, Y. Jiang nd Z.-H. Zhou. Partial Multi-View Clustering. In: Proceedings of the 28th AAAI Conference on
% Artificial Intelligence (AAAI'14),Qu��bec, Canada ,2014.

% Contact handong.zhao@gmail.com if you have any questions
% Zhao, et al., Incomplete Multi-modal Visual Data Grouping, IJCAI'16

%%
if (min(truth)==0)
    truth = truth + 1;
end

option.option = numClust;
option.truth = truth;
[Ux,Uy,P2,P1,P3,W] = IMG(X2,Y2,X1,Y3,option);

% fprintf('running spectral clustering...\n');
kmeans_avg_iter = 10;
for i=1: kmeans_avg_iter

    C = clu_ncut(W,numClust);
    C = C';
    
    %%
    [A nmii(i) avgenti(i)] = compute_nmi(truth,C);
    [Fi(i),Pi(i),Ri(i)] = compute_f(truth,C);
    [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C);
end
F = mean(Fi);
P = mean(Pi);
R = mean(Ri);
nmi = mean(nmii);
avgent = mean(avgenti);
AR = mean(ARi);

fprintf('nmi: %f(%f)\n', nmi, std(nmii));

