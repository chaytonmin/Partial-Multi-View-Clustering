function jacc = Jac(K,knum,KF)

n = size(K,1);
jacc = 0;
for i = 1:n
    setI = setdiff([1:n],i);
    KFi = KF(i,setI);
    [val1,indx1] = sort(KFi,'descend');
    nbikf = indx1(1:knum);
    
    Ki = K(i,setI);
    [val2,indx2] = sort(Ki,'descend');
    nbik = indx2(1:knum);
    
    nbintersec = length(intersect(nbikf,nbik));
    nbunion   = length(union(nbikf,nbik));
    
    jacc = jacc + (1/n)*(nbintersec/nbunion);
end

