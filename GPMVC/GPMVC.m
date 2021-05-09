function [U,centroidPc, P,objValue] = GPMVC(Xpaired, Xsingle,W,views,options)
%      partial view multi-view data set :   X2: n2*dx  Y2: n2*dy  X1: n1*dx  Y3: n3*dy
%      (X2,Y2); examples appearing in both views; 
%      (X1,);   examples appearing in only view x
%      (,Y3):   examples appearing in only view y  
%      W1 : Weight matrix for the affinity graph of View 1
%      W2 : Weight matrix for the affinity graph of View 2
%      dx,dy:   feature dimension of viewx viewy
%      option.alpha: parameter controls the importance of views
%      option.latentdim: the feature dimension of latent space(default value: cluster number)
%   ouptut: 
%     Ux, Uy:   (k*dx, k*dy)  basis of latent space
%     P2,P1,P3: (n2*k, n1*k, n3*k)  data representations for examples in the latent space 
%     centroidP2 : Refers to the consensus matrix for the complete datapoints
%	Writen by Jialu Liu (jliu64@illinois.edu)
% 	Modified by Zhenfan Wang (zfwang@mail.dlut.edu.cn)
%   Further modified by Nishant Rai (nishantr AT iitk DOT ac DOT in)

%Note that columns are data vectors here
%Number of views are assumed to be 2 in this code

tic;
Rounds = options.rounds;
alpha=options.alpha;
beta=options.beta;
error=options.error;
gamma = options.gamma;

K = options.K;

A=cell(views,1);
L=cell(views,1);
numInst=zeros(views,1);
Feat=zeros(views,1);
singX=zeros(views,1);
for v=1:views
    A{v}=horzcat(Xpaired{v}',Xsingle{v}');
    [numInst(v), Feat(v)]=size(A{v}');
    P{v}=rand(numInst(v),K); 
    U=rand(Feat(v),K);
    nSmp = size(W{v},1);
    Wtemp = W{v};  
    DCol = full(sum(Wtemp,2));
    D = spdiags(DCol,0,nSmp,nSmp);
    L{v} = D - Wtemp;  
    singX(v)=size(Xsingle{v}',1);
end  
%A1 = horzcat(X1,X2);
%A2 = horzcat(Y2,Y3);

%[numInst1,Featx]=size(A1');                             %Number of instances, FeatX: Number of features in view 1
%[numInst2,Featy]=size(A2');                             %Number of instances with view 2

%P1=rand(numInst1,K);                                %Random initialization
%P2=rand(numInst2,K);
%Ux=rand(Featx,K);
%Uy=rand(Featy,K);

%nSmp = size(W(1),1);
%Wtemp = W1;            %Modify the weight matrix with the involved parameters
%DCol = full(sum(Wtemp,2));
%D = spdiags(DCol,0,nSmp,nSmp);
%L1 = D - Wtemp;                              %Get matrix L

%nSmp = size(W2,1);
%Wtemp = W2;            %Modify the weight matrix with the involved parameters
%DCol = full(sum(Wtemp,2));
%D = spdiags(DCol,0,nSmp,nSmp);
%L2 = D - Wtemp;                              %Get matrix L

%singX=size(X1',1);
numCom=size(Xpaired{1},1);

objValue = [];

weights = ones(1,views)./views;

%% initialize basis and coefficient matrices, initialize on the basis of standard GNMF algorithm
tic;
Goptions.alpha = options.alpha*options.beta;
rand('twister',5489);
U=cell(views,1);
P=cell(views,1);
for v=1:views
[U{v}, P{v}] = GNMF(A{v}, K, W{v}, Goptions); 
[U{v}, P{v}] = Normalize(U{v}, P{v});
end%In this case, random inits take place
%rand('twister',5489);
%[Uy, P2] = GNMF(A2, K, W2, Goptions);
toc;
%%

%[Ux, P1] = Normalize(Ux, P1);
%[Uy, P2] = Normalize(Uy, P2);

%% Alternate Optimisations for consensus matrix and individual view matrices
optionsPGNMF = options;
j = 0;
sumRound=0;
H=cell(views,1);
while j < Rounds                            %Number of rounds of AO
    sumRound=sumRound+1;
    j = j + 1;
    
    %centroidPc = (weights(1)^gamma)*(options.alphas(1)*P1(singX+1:end,:)) + (weights(2)^gamma)*(options.alphas(2)*P2(1:numCom,:));     
     centroidPc=0;
     tmpsum=0;
   for v=1:views
         centroidPc=  centroidPc+(weights(v)^gamma)*(options.alphas(v)*P{v}(1:numCom,:));
         tmpsum=tmpsum+(weights(v)^gamma)*options.alphas(v) ;
   end
    
    %From the paper, we have a definite solution for Pc*
  %  tmpSum = (weights(1)^gamma)*options.alphas(1) + (weights(2)^gamma)*(options.alphas(2));
    centroidPc = centroidPc / tmpsum;
        
    %Update the weights if the corresponding option is set
    if (options.varWeight > 0)
        
        tmpSum=0;
        for  v=1:views
        tmp1 = (A{v} - U{v}*P{v}');
        tmp2 = (P{v}(singX(v)+1:end,:) - centroidPc);
        H{v}= gamma*(sum(sum(tmp1.^2)) + options.alphas(v)*(sum(sum(tmp2.^2)))+ (beta*alpha)*sum(sum((P{v}'*L{v}).*P{v}')));
         H{v} = H{v}^(1.0/(1-gamma));
          tmpSum = tmpSum+H{v};
        end
       for v=1:views
           weights(v)=(H{v}/tmpSum);
       end
       
    end
    
   % for i=1:views
      %  fprintf('%.3f ',weights(i));
   % end
    %fprintf(' weights: %d\n',j);
    
    logL = 0;                                   %Loss for the round
    
    %Compute the losses
    
    for v=1:views
       tmp1 = (A{v} - U{v}*P{v}');
       tmp2 = (P{v}(1:numCom,:) - centroidPc);
        logL = logL + sum(sum(tmp1.^2)) + options.alphas(v)*(sum(sum(tmp2.^2)))+ (beta*alpha)*sum(sum((P{v}'*L{v}).*P{v}'));
    end
    
 
   % fprintf('LOG %.9f, ',logL);
    objValue = [objValue logL];                %End indicates last index of array, so basically push operation
    
   % if mod(j,10)==0
   % fprintf('Iteration %d, objective value %g\n', j, objValue(j));
    %end

    if j>1 && ((abs(objValue(j)-objValue(j-1))/objValue(j) < error)|| objValue(j)<=error)
        fprintf('Objective value converge to %g at iteration %d before the maxIteration reached \n',objValue(j),j);
        break;
    end
    
    if sumRound==Rounds
        break;
    end

    %Peform optimization with Pc* (centroidPc) fixed and inits finalU, finalV
    
   optionsPGNMF.begins = 1;
    optionsPGNMF.ends = numCom; 
    for v=1:views
    optionsPGNMF.alphaPriv = options.alphas(v);
    Ptmp = [centroidPc;P{v}(numCom+1:end,:)];
    [U{v}, P{v}] = PartialGNMF(A{v}, K, Ptmp, W{v}, optionsPGNMF, U{v}, P{v});
    %W has not been multiplied by the weight
    end
    

end

%tmp2 = (P1(singX+1:end,:) - centroidPc);
%sum(sum(tmp2.^2))
%tmp2 = (P2(1:numCom,:) - centroidPc);
%sum(sum(tmp2.^2))   


%P1 = P1(1:singX,:);
%P2 = P2(numCom+1:end,:);
for v=1:views
P{v}= P{v}(numCom+1:end,:);
end

fprintf('\n');

%workspace
toc


function [U, V] = Normalize(U, V)
    [U,V] = NormalizeUV(U, V, 0, 1);

function [U, V] = NormalizeUV(U, V, NormV, Norm)
    nSmp = size(V,1);
    mFea = size(U,1);
    if Norm == 2
        if NormV
            norms = sqrt(sum(V.^2,1));
            norms = max(norms,1e-10);
            V = V./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);
        else
            norms = sqrt(sum(U.^2,1));
            norms = max(norms,1e-10);
            U = U./repmat(norms,mFea,1);
            V = V.*repmat(norms,nSmp,1);
        end
    else
        if NormV
            norms = sum(abs(V),1);
            norms = max(norms,1e-10);
            V = V./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);
        else
            norms = sum(abs(U),1);
            norms = max(norms,1e-10);
            U = U./repmat(norms,mFea,1);
            V = bsxfun(@times, V, norms);
        end
    end
