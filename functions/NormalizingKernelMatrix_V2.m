function NK = NormalizingKernelMatrix_V2(K)
Kd = diag(diag(K) + 10^(-6));
% Zero_Index = Kd < 10^(-6);
% Kd(Zero_Index) = 1;
% CKd = repmat(Kd , [1 , length(Kd)]);
% RKd = repmat(Kd' , [length(Kd) , 1]);
% K = K./CKd;
% K = K./RKd;
NK = K / Kd;
end