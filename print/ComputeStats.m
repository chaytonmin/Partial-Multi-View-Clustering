function [result] = ComputeStats(X, label, K, kmeansIter, kmeansFlag)
    if (~exist('kmeansFlag','var'))
        kmeansFlag = 1;
    end
    if (~exist('kmeansIter','var'))
        kmeansIter = 10;
    end
    if (kmeansFlag == 1)
%     		fprintf('running k-means... ');
    end
    for i=1:kmeansIter
        rand('twister',5489);
        if kmeansFlag > 0
            indic = litekmeans(X, K, 'Replicates',20);
        else
            [~, indic] = max(X, [] ,2);
        end
        index = bestMap(label, indic);
        
        [G]=cmat(label, index);
        imagesc(G);
        
        [ac(i), nmi_value(i), cnt(i)] = CalcMetrics(label, indic);
        [Pri(i)] = purity(label, indic);
        [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(label,indic);
        [Fi(i),Pi(i),Ri(i)] = compute_f(label,indic);
    end
    if (kmeansFlag == 1)
%         fprintf('ac: %0.4f,  nmi:%0.4f,  purity:%.4f,  errors: %d/%d\n', ...
%                         mean(ac), mean(nmi_value), mean(Pri), round(mean(cnt)), length(label));    
    end
    ac = mean(ac);
    nmi=mean(nmi_value);
    F = mean(Fi);
    P = mean(Pi);
    R = mean(Ri);
    AR = mean(ARi);
    Pur=mean(Pri);
    result=[ac,nmi,AR,F,P,R,Pur];
    %avgent = mean(avgenti);
    disp(sprintf('ACC: %0.4f  NMI:%0.4f F-score:%0.4f Precision:%0.4f Recall:%0.4f AR:%0.4f Purity:%0.4f \n', ac,nmi,F,P,R,AR,Pur));
end
