clear
clc
warning off;
addpath(genpath('ClusteringEvaluation/'));
addpath(genpath('lmkkmeans-master/'));

path = '';
addpath(genpath(path));
dataName = 'yale-0.1-Kernel-MKKIK';
%% %% washington; wisconsin; texas; cornell
load([path,'datasets/',dataName,'.mat'],'KH','Y','S');%KH：核，Y：类标,S：缺失的哪个数据点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qnorm = 2;
% [H_normalized,gamma,obj] = mkkmeans_train(KH,numclass,qnorm);
% res_gnd = myNMIACC(H_normalized,Y,numclass);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsionset =[0.9];
iter = 1;
for ie =1:length(epsionset)
   
        tic;
        algorithm_choose3 = 'algorithm0';
        [H_normalized,gamma,obj,KH] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose3);
        timingcost = toc;
        res = myNMIACC(H_normalized,Y,numclass);
         fprintf('Acc: %.4f  NMI: %.4f AR: %.4f F: %.4f  P: %.4f R: %.4f Purity: %.4f\n',res); 
        save([path,'work2016/myFinalRes/',dataName,'_missingRatio_',num2str(epsionset),'_norm_',num2str(qnorm),...
            '_clustering_iter_',num2str(iter),'.mat'],'res','timingcost');
      
end