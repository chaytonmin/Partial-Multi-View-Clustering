function ka=algorithm0(K,S)
%È±Ê§µÄºËÌî³ä0
num = size(K,1);
nbkernel = size(K,3);
ka = zeros(num,num,nbkernel);
 for p =1:nbkernel
      mis_indx = S{p}.index;
      obs_indx = setdiff(1:num,mis_indx);  
      ka(obs_indx,obs_indx,p)=K(obs_indx,obs_indx,p);
 end
