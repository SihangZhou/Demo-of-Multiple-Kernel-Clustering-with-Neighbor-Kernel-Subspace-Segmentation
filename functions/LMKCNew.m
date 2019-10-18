function [Mu , Z , flag , TotalObj] = LMKCNew(KHL , M1 ,  Alpha_C , Beta_C , Gamma_C , L , LargestIteration)

[SampleNum , ~ , KerNum] = size(KHL);

Mu = ones(KerNum , 1) / KerNum;
TotalObj = zeros(LargestIteration , 1);
Iter = 0;
flag = 0;
while Iter < LargestIteration && flag == 0
    CurrentKernel = sumKbeta(KHL , Mu);
    % Calculate Z
    KSquare = CurrentKernel*CurrentKernel;
    A = KSquare + Alpha_C*eye(SampleNum) + Beta_C * L;
    AInv = eye(size(A)) / A;
    %     AInv = inv(A);

    Z = AInv * KSquare;
    
    % Calcluate Mu
    M2 = M_2Calculation_New(KHL , Z);
    Mu = MuCulationNew(M1 , M2 , Gamma_C);
    
    Obj = ObjCalculation(Mu , M1 , M2 , Alpha_C, Beta_C, Gamma_C , Z , L);
    Iter = Iter + 1;
    TotalObj(Iter) = Obj;
    
    if Iter >=2
        if (abs( TotalObj(Iter)-TotalObj(Iter-1) ) < TotalObj(Iter) * 10^(-5))
            flag = 1;
        end
    end
end
TotalObj  = TotalObj(1 : Iter);