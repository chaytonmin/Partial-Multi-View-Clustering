clear 
clc
addpath(genpath('measure/'));
addpath(genpath('misc/'));
addpath(genpath('print/'));
datasetdir='data/';
resultdir='Results/';
load L0.01.mat;
%load U0.01.mat;
%load Yc0.01.mat;
dataname={'animals'};ve='L';
per=0.01;
maxiter=100;
para=[0.001,0.003,0.005,0.007,0.009,0.01];
k=50;
%learning 
for id=1:length(per)
    dataf=strcat(datasetdir,dataname,num2str(per(id)),'.mat');
    datafname=cell2mat(dataf(1));
    load (datafname);
    views=length(Xpaired);
    npaired=size(Xpaired{1,1},1);
    Yc0=rand(npaired,k);
    
 % compute L
   for i=1:views
       Xpaired{i}=Xpaired{i}';
       Xsingle{i}=Xsingle{i}';
        Xa{i}=[Xpaired{i},Xsingle{i}];
        gnd_all{i}=[truthPaired;truthSingle{i}];
        [dims(i),nsingle(i)]=size(Xsingle{i});
        Y0{i}=rand(nsingle(i),k);
        U0{i}=rand(dims(i),k)
   end
     gamma=para(4);
     beta=para(3);
    L=computeL(Xa,gnd_all);
   tic
    [U,Yall,obj, acc_sym] =unified_subspace(L, Xpaired, Xsingle, Xa,gamma, beta,maxiter, Yc0, Y0,U0,truthF);
    time=toc;
    fmodel2=sprintf('Resultes/%s.mat',dataset);
save(fmodel2,'acc_sym','Time','U','Yall');
end