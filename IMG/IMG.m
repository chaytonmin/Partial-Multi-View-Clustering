function [U1,U2,P2,P1,P3,A] = IMG(X2,Y2,X1,Y3,option)
% by handong Zhao
% contact handong.zhao@gmail.com if you have any questions
% Zhao, et al., Incomplete Multi-modal Visual Data Grouping, IJCAI'16

warning off;
lamda = option.lamda;
beta = option.beta;
gamma = option.gamma;
k = option.latentdim;
numClust = option.option;
truth = option.truth;

[numInst1,Featx]=size(X1);
[numInst2,Featx]=size(X2);
[numInst3,Featy]=size(Y3);

P1init=rand(numInst1,k);
P2init=rand(numInst2,k);
P3init=rand(numInst3,k);
Uxinit=rand(k,Featx);
Uyinit=rand(k,Featy);

P1=P1init;
P2=P2init;
P3=P3init;

Q1=P1init;
Q2=P2init;
Q3=P3init;

U1=Uxinit;
U2=Uyinit;


tol = 1e-8;
maxIter = 200;
mu = 1e-1;
rho = 1.1;
mu_bar = 1e10;


% initialize Y
L1 = zeros(size(P1));
L2 = zeros(size(P2));
L3 = zeros(size(P3));

I = eye(k);
numInst = numInst1+ numInst2+ numInst3;
OneN = ones(numInst,1);
In = eye(numInst);
La = rand(numInst,numInst);
oldA = zeros(size(numInst,numInst));

for ii = 1: maxIter
    
    % update U1  ------------------------
    U1 = (P2'*P2+lamda*I)\(P2'*X2);
    

    % update U2  ------------------------
    U2 = (P2'*P2+lamda*I)\(P2'*Y2);
    
    % update P2
    tmpA = 2*X2*U1'+2*Y2*U2'-L2+mu*Q2;
    tmpB = 2*U1*U1'+2*U2*U2'+mu*I;
    P2 = tmpA/(tmpB);
    
    
    % update P1  ------------------------
    tmpA = 2*X1*U1'-L1+mu*Q1;
    tmpB = 2*U1*U1'+mu*I;
    P1 = tmpA/(tmpB);
    
    % update P3  ------------------------
    tmpA = 2*Y3*U2'-L3+mu*Q3;
    tmpB = 2*U2*U2'+mu*I;
    P3 = tmpA/(tmpB);
    
    
    
    % update Q = [Q2; Q1; Q3]
    L = [L2; L1; L3];
    P = [P2; P1; P3];
    Q = pinv(beta*(La'+La)+mu*In)*(L+mu*P);
    Q2 = Q(1:numInst2,:);
    Q1 = Q(numInst2+1:numInst1+numInst2,:);
    Q3 = Q(numInst1+numInst2+1:end,:);
    
    % update A
    
    for iSamp = 1:size(Q,1),
        tmpVal = 0;
        Qi = Q(iSamp,:);
        di = zeros(numInst,1);
        for jSamp = 1:size(Q,1),
            Qj = Q(jSamp,:);
            tmpVal = tmpVal + beta/(4*gamma)*norm((Qi-Qj),'fro')^2;
            di(jSamp) = beta/(4*gamma)*norm((Qi-Qj),'fro')^2;
        end
        A(:,iSamp) = ((1+ tmpVal)/numInst)*OneN - di;
        A = max(A,0);
    end
    
    A = (A+A')/2;
    
    if ii==1,
        resA = zeros(size(A));
    else
        resA = (A - oldA);
    end
    oldA = A;
    
    La = diag(sum(A,2))- A;
    
    leq1 = P1-Q1;
    leq2 = P2-Q2;
    leq3 = P3-Q3;
    L1 = L1 + mu * leq1;
    L2 = L2 + mu * leq2;
    L3 = L3 + mu * leq3;
    mu = min(mu*rho, mu_bar);
    
    % convergence criterion
    stopCriterion = max(max(max(abs(leq1))),max(max(abs(leq2))));
    stopCriterion = max(stopCriterion,max(max(abs(leq3))));
    
    if mod(ii, 10) == 0
        C = clu_ncut(A,numClust);
        C = C';
        [~, nmi_sc, ~] = compute_nmi(truth,C);
               
        UPI=[P2;P1;P3];
        if (1)
            norm_mat = repmat(sqrt(sum(UPI.*UPI,2)),1,size(UPI,2));
            %%avoid divide by zero
            for i=1:size(norm_mat,1)
                if (norm_mat(i,1)==0)
                    norm_mat(i,:) = 1;
                end
            end
            PIn = UPI./norm_mat;
        end
        for i=1: 20
            C = kmeans(PIn,numClust,'EmptyAction','drop');
            [~, nmii(i), ~] = compute_nmi(truth,C);
        end
        nmi_kmeans = mean(nmii);
        
        
        
        disp(['iter ' num2str(ii) ', para ' num2str(lamda) '/' num2str(beta) '/' num2str(gamma)...
            ', stopC ' num2str(stopCriterion) ', NMI_SC ' num2str(nmi_sc)...
            ', NMI_kmeans ' num2str(nmi_kmeans)]);
    end
    
    if stopCriterion < tol
        break;
    end
    
end
end