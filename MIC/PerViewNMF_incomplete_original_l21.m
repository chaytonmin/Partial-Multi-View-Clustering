function [U_final, V_final, nIter_final, elapse_final, bSuccess, objhistory_final] = PerViewNMF_incomplete_original_l21(X, k, Uo, options, U, V, C)
%
% Notation:
% X ... (nSmp x mFea) data matrix of one view
%       mFea  ... number of features
%       nSmp  ... number of samples
% k ... number of hidden factors
% Uo... consunsus
% options ... Structure holding all settings
% U ... initialization for coefficient matrix
% V ... initialization for basis matrix
%
%   Originally written by Deng Cai (dengcai AT gmail.com) for GNMF
%   Modified by Weixiang Shao (wshao4@uic.edu)

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIterOrig = options.minIter;
minIter = minIterOrig-1;
meanFitRatio = options.meanFitRatio;

alpha = options.alpha;
beta = options.beta;

[nSmp, mFea]=size(X);

bSuccess.bSuccess = 1;

selectInit = 1;
if isempty(U)
    U = abs(rand(nSmp,k));
    V = abs(rand(mFea,k));
else
    nRepeat = 1;
end

[U,V] = Normalize(U, V);
if nRepeat == 1
    selectInit = 0;
    minIterOrig = 0;
    minIter = 0;
    if isempty(maxIter)
        objhistory = CalculateObj(X, U, V, Uo, C, alpha, beta);
        meanFit = objhistory*10;
    else
        if isfield(options,'Converge') && options.Converge
            objhistory = CalculateObj(X, U, V, Uo, C, alpha, beta);
        end
    end
else
    if isfield(options,'Converge') && options.Converge
        error('Not implemented!');
    end
end



tryNo = 0;
while tryNo < nRepeat
    tmp_T = cputime;
    tryNo = tryNo+1;
    nIter = 0;
    maxErr = 1;
    nStepTrial = 0;
    %disp a
    while(maxErr > differror)
        % ===================== update U ========================

        XV = (C.^2)*X*V;  % mnk or pk (p<<mn)
        VV = V'*V;  % mk^2
        UVV = (C.^2)*U*VV; % nk^2

        XV = XV + alpha * (C.^2) * Uo;
        D = zeros(size(U,1));
        for j = 1:size(U,1)
            D(j,j) = 1/norm(U(j,:));
        end
        UVV = UVV + alpha * (C.^2) * U + 0.5*beta*D*U;

        U = U.*sqrt(XV./max(UVV,1e-30));

        % ===================== update V ========================

        XU = X'*(C.^2)*U;
        UU = U'*(C.^2)*U;
        VUU = V*UU;
        V = V.*sqrt(XU./max(VUU,1e-30));

        [U,V] = Normalize(U, V);
        nIter = nIter + 1;
        if nIter > minIter
            if selectInit
                objhistory = CalculateObj(X, U, V, Uo, C, alpha, beta);
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj = CalculateObj(X, U, V, Uo, C, alpha, beta);
                    objhistory = [objhistory newobj];
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        newobj = CalculateObj(X, U, V, Uo, C, alpha, beta);
                        objhistory = [objhistory newobj];
                    end
                    maxErr = 1;
                    if nIter >= maxIter
                        maxErr = 0;
                        if isfield(options,'Converge') && options.Converge
                        else
                            objhistory = 0;
                        end
                    end
                end
            end
        end
    end

    elapse = cputime - tmp_T;

    if tryNo == 1
        U_final = U;
        V_final = V;
        nIter_final = nIter;
        elapse_final = elapse;
        objhistory_final = objhistory;
        bSuccess.nStepTrial = nStepTrial;
    else
       if objhistory(end) < objhistory_final(end)
           U_final = U;
           V_final = V;
           nIter_final = nIter;
           objhistory_final = objhistory;
           bSuccess.nStepTrial = nStepTrial;
           if selectInit
               elapse_final = elapse;
           else
               elapse_final = elapse_final+elapse;
           end
       end
    end

    if selectInit
        if tryNo < nRepeat
            %re-start
            U = abs(rand(mFea,k));
            V = abs(rand(nSmp,k));
            [U,V] = Normalize(U, V);
        else
            tryNo = tryNo - 1;
            minIter = 0;
            selectInit = 0;
            U = U_final;
            V = V_final;
            objhistory = objhistory_final;
            meanFit = objhistory*10;

        end
    end
end

nIter_final = nIter_final + minIterOrig;
[U_final, V_final] = Normalize(U_final, V_final);
end

%==========================================================================

function [obj, dV] = CalculateObj(X, U, V, Uo, C, alpha, beta)
    tmp = C*(U-Uo);
    obj_Lap = sum(sum(tmp.^2));
    dX = C*(U*V'-X);
    obj_NMF = sum(sum(dX.^2));
    obj_L1 = sum(sum(abs(U)));
    tmp3 = 0;
        for k =1:size(U,2);
            tmp3 = tmp3 + norm(U(:,k));
        end
    obj = obj_NMF+ alpha * obj_Lap + beta * tmp3;
end

function [U, V] = Normalize(U, V)
    nSmp = size(U,1);
    mFea = size(V,1);
	norms = sum(abs(V),1);
    norms = max(norms,1e-30);
    V = V./repmat(norms,mFea,1);
    U = bsxfun(@times, U, norms);
end
