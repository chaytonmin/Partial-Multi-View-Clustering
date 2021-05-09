function indx_0 = genarateNeighborhood(KC,tau)

num = size(KC,1);
KC0 = KC - 10^8*eye(num);
[val,indx] = sort(KC0,'descend');
indx_0 = indx(1:tau,:);