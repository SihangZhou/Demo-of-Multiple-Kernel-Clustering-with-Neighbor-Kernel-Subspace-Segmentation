function M = calculateM(K)

numker = size(K,3);
M = zeros(numker);
for p =1:numker
    for q = p:numker
        M(p,q) = trace(K(:,:,p)'*K(:,:,q));
    end
end
M = (M+M')-diag(diag(M));