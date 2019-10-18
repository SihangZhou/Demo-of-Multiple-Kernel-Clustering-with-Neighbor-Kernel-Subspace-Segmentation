function M_2 = M_2Calculation_New(Ktrntrn , Current_A)

ZI = Current_A - eye(size(Current_A));  % Calculate Z-I

TempKtrntrn = zeros(size(Ktrntrn));
KernelNum = size(Ktrntrn , 3);

for IKernel = 1 : KernelNum
    TempKtrntrn(: ,: , IKernel) = Ktrntrn(: ,: , IKernel) * ZI;
end
clear Ktrntrn

M_2 = zeros(KernelNum, KernelNum);
for KNum = 1 : KernelNum
    for KNum2 = KNum : KernelNum
        M_2(KNum , KNum2) = trace(TempKtrntrn(: , : , KNum) * TempKtrntrn(:, : , KNum2)');
        M_2(KNum2 , KNum) = M_2(KNum , KNum2);
    end
end

end