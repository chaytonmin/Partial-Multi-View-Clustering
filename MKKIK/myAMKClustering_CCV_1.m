clear
clc
warning off;
addpath(genpath('ClusteringEvaluation/'));
addpath(genpath('lmkkmeans-master/'));

path = '';
addpath(genpath(path));
dataName = 'washington'; %%% flower17; flower102; CCV; caltech101_numofbasekernel_10
%% %% washington; wisconsin; texas; cornell
load([path,'datasets/',dataName,'_Kmatrix'],'KH','Y');
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
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
epsionset = [0.1:0.1:0.9];
for ie =1:length(epsionset)
     for iter = 1:10 
        load([path,'work2016/generateAbsentMatrix/',dataName,'_missingRatio_',num2str(epsionset(ie)),...
            '_missingIndex_iter_',num2str(iter),'.mat'],'S');
        
%         %%%%%%%%%%%--Zero-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
%         KH1 = algorithm2(KH,S);
%         [H_normalized1,gamma1,obj1] = mkkmeans_train(KH1,numclass,qnorm);
%         timingcost(1) = toc;
%         res(:,1) = myNMIACC(H_normalized1,Y,numclass);
%         
% %         %%%%%%%%%%%--mean-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         KH2 = algorithm3(KH,S);
%         [H_normalized2,gamma2,obj2] = mkkmeans_train(KH2,numclass,qnorm);
%         timingcost(2) = toc;
%         res(:,2) = myNMIACC(H_normalized2,Y,numclass);
%         
%         
%         %%%%%%%%%%--knn-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         KH3 = algorithm0(KH,S,7);
%         [H_normalized3,gamma3,obj3] = mkkmeans_train(KH3,numclass,qnorm);
%         timingcost(3) = toc;
%         res(:,3) = myNMIACC(H_normalized3,Y,numclass);
%         
%         % %         %%%%%%%%%%%---EM-filling---%%%%%%%%%%%%%%%%%%%%%%%
%         %         KH4 = algorithm6(KH,S);
%         %         [H_normalized4,gamma4,obj4] = mkkmeans_train(KH4,numclass,qnorm);
%         %         res(:,4) = myNMIACC(H_normalized4,Y,numclass);
%         %%%%%%%%%--Laplacian-filling----%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         alpha04 = 1e-3;
%         KH4 = algorithm4(KH,S,numclass,alpha04);
%         [H_normalized4,gamma4,obj4] = mkkmeans_train(KH4,numclass,qnorm);
%         timingcost(4) = toc;
%         res(:,4) = myNMIACC(H_normalized4,Y,numclass);
% %         
% %         %%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         algorithm_choose1 = 'algorithm2';
%         [H_normalized5,gamma5,obj5,KH5] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose1);
%         timingcost(5) = toc;
%         res(:,5) = myNMIACC(H_normalized5,Y,numclass);
%         
%         tic;
%         algorithm_choose2 = 'algorithm3';
%         [H_normalized6,gamma6,obj6,KH6] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose2);
%         timingcost(6) = toc;
%         res(:,6) = myNMIACC(H_normalized6,Y,numclass);
%         
%         tic;
        algorithm_choose3 = 'algorithm0';
        [H_normalized7,gamma7,obj7,KH7] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose3);
        timingcost(7) = toc;
        res(:,7) = myNMIACC(H_normalized7,Y,numclass);
        

        
% % %         %%%%%%%%---IJCAI2017-----%%%%%%%%%%%%%%%%
%         algorithm_choose8 = 'algorithm3';
%         lambdaset8 = 2.^[-15:2:11];
%         accval8 = zeros(length(lambdaset8),1);+-
%         nmival8 = zeros(length(lambdaset8),1);
%         purval8 = zeros(length(lambdaset8),1);
%         algval8 = zeros(length(lambdaset8),1);
%         for il =1:length(lambdaset8)
%             tic;
%             [H_normalized8,gamma8,obj8,KH8] = myamkcwithlambda(KH,S,numclass,qnorm,algorithm_choose8,lambdaset8(il));
%             timingcost(8) = toc;
%             res8 = myNMIACC(H_normalized8,Y,numclass);
%             accval8(il) = res8(1);
%             nmival8(il) = res8(2);
%             purval8(il) = res8(3);
%             algval8(il) = calKernelAlignment(KH,KH8)'*gamma8;
%         end
%         res(:,8) = [max(accval8); max(nmival8); max(purval8)];
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         alignment(1) = calKernelAlignment(KH,KH1)'*gamma1;
%         alignment(2) = calKernelAlignment(KH,KH2)'*gamma2;
%         alignment(3) = calKernelAlignment(KH,KH3)'*gamma3;
%         alignment(4) = calKernelAlignment(KH,KH4)'*gamma4;
%         alignment(5) = calKernelAlignment(KH,KH5)'*gamma5;
%         alignment(6) = calKernelAlignment(KH,KH6)'*gamma6;
%         alignment(7) = calKernelAlignment(KH,KH7)'*gamma7;
%         alignment(8) = max(algval8);
%         
        save([path,'work2016/myFinalRes/',dataName,'_missingRatio_',num2str(epsionset(ie)),'_norm_',num2str(qnorm),...
            '_clustering_iter_',num2str(iter),'.mat'],'res','timingcost','alignment');
      end
end