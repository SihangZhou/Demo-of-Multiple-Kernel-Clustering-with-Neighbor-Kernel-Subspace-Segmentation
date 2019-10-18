function [KHL] = LocalKernelCalculation(KH  , NNRate)

[SampleNum , ~ , KerNum] = size(KH);

gamma0 = ones(KerNum,1)/KerNum;
avgKer  = mycombFun(KH,gamma0);

KHL = zeros(size(KH));
for kk = 1 : KerNum
    KHL(:,:,kk) = NeighborKernelSH(KH(:,:,kk) , avgKer, NNRate);
    KHL(:,:,kk) = ( KHL(:,:,kk) + KHL(:,:,kk)' ) / 2;
end

for BKNum = 1 : KerNum
    KHL(:,:,BKNum) = NormalizingKernelMatrix_V2(KHL(:,:,BKNum));
    KHL(:,:,BKNum) = (KHL(:,:,BKNum) + KHL(:,:,BKNum)')/2 + eye(SampleNum)*10^(-8);
end
end