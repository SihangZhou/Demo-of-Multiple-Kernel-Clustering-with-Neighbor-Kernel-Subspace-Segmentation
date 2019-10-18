% function [M , M2] = MCalculation_New(K)
function [M] = MCalculation_New(K)
% This function calculate the matrix M in paper
KernelNum = size(K , 3);
M = zeros(KernelNum, KernelNum);
% M2 = zeros(KernelNum, KernelNum);
for KNum = 1 : KernelNum
    for KNum2 = KNum : KernelNum
        M(KNum , KNum2) = trace(K(: , : , KNum)' * K(:, : , KNum2));
        M(KNum2 , KNum) = M(KNum , KNum2);
    end
end

% for KNum = 1 : KernelNum
%     for KNum2 = KNum : KernelNum
%         M2(KNum , KNum2) = M(KNum , KNum2) / sqrt(M(KNum2 , KNum2) * M(KNum , KNum));
%         M2(KNum2 , KNum) = M2(KNum , KNum2);
%     end
% end


end
