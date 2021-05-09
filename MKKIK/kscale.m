function K = kscale(K)

for i=1:size(K,3)
    K(:,:,i) = K(:,:,i) ./ sqrt(diag(K(:,:,i))*diag(K(:,:,i))');
end
