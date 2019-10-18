clear
clc
warning off;

DataName = cell(10 , 1);
DataName{1} = 'YALE';
DataName{2} = 'caltech101_numOfClass_10';
DataName{3} = 'proteinFold';
DataName{4} = 'UCI_DIGIT';
DataName{5} = 'flower17';
DataName{6} = 'psortPos';
DataName{7} = 'flower102_numOfClass_90';
DataName{8} = 'caltech101_10base';
DataName{9} = 'plant';
DataName{10} = 'flower102_numOfClass_100';
DataName{11} = 'caltech101_numOfClass_40';
DataName{12} = 'Reuters';
DataName{13} = 'psortNeg';


for ICount = [1 : 7 , 9 , 10 , 11 , 12]
    
    res = zeros(3 , 8);
    
    path = 'D:\Others\sihangzhou_code\18_GlobalLocalMKC\1718Final_Combination\';
    pathdata = 'D:\Others\sihangzhou_code\17_SpectralKMeans\datasets_xinwang16';
    addpath(genpath(path));
    dataName = DataName{ICount}; %%% flower17; flower102; CCV; caltech101_mit_numOfClass_10
    %% caltech101_numOfClass_10_Kmatrix
    %% washington texas cornell; wisconsin, UCI_DIGIT
    load([pathdata,'\',dataName,'_Kmatrix'],'KH','Y');
    % load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CluNum = length(unique(Y));
    KerNum = size(KH,3);
    SampleNum = size(KH,1);
    LowRankNum = round(5 * log2(SampleNum));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KH = kcenter(KH);
    KH = knorm(KH);
    M = calculateM(KH);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qnorm = 2;
    %%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma0 = ones(KerNum,1)/KerNum;
    avgKer  = mycombFun(KH,gamma0);
    [H_normalized1] = mykernelkmeans(avgKer,CluNum);
    res(:,1) = myNMIACC(H_normalized1,Y,CluNum);
    
    % %%%%%%%%%%---Single Best%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    accpval2 = zeros(KerNum,1);
    nmipval2 = zeros(KerNum,1);
    purpval2 = zeros(KerNum,1);
    for p =1:KerNum
        [H_normalized2] = mykernelkmeans(KH(:,:,p),CluNum);
        res2 = myNMIACC(H_normalized2,Y,CluNum);
        accpval2(p) = res2(1);
        nmipval2(p) = res2(2);
        purpval2(p) = res2(3);
    end
    res(:,2) = [max(accpval2);max(nmipval2);max(purpval2)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%       GLMKC - without LLE      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  algorithm set up && super parameter setting
    [KHL] = LocalKernelCalculationNew(KH , LowRankNum);
    LargestIteration = 15;
    KernelWeight = ones(1 , KerNum) / KerNum;
    avgKer=sumKbeta(KHL,KernelWeight);
    BaseValue = norm(avgKer , 'fro')^2;
    RegularizationValue = 10^(-4) * BaseValue;
    Alpha_T = 2 * 0;
    Beta_T = 2.^(-4 : 4 : 4);
    NNRate = [0.01];
    
    GLMKC_ResZ_NLLE = zeros(length(NNRate) , length(Beta_T) , 3);
    GLMKC_ResK_NLLE = zeros(length(NNRate) , length(Beta_T) , 3);
    
    %% Algorithm
    for IRate = 1 : length(NNRate)
        
        % Preprocessing (Calculating KHL, L , M1)
        CRate = NNRate(IRate);
        [KHL , L , M1] = LGMKCPreprocessFinalNew(KH , LowRankNum);
        
        for IBeta = 1 : length(Beta_T)
            Beta_C = Beta_T(IBeta);
            
            InputKnum = max(2* CluNum , LowRankNum);
            
            % Kernel weight and Z calculation
            [Mu , Z , flag , TotalObj] = LGMKC(KHL, M1, RegularizationValue , Alpha_T , Beta_C , L , LargestIteration , InputKnum);
            
            PON = TotalObj(1:end-1) - TotalObj(2 : end);
            DON = sum(PON < -10^(-10)) == 0;
            fprintf('\n\n The Current dataset is : %d; IBeta = %d , IGamma = %d , flag = %d \n\n ' , ICount , IRate , IBeta , DON)
            
            % Result Record
            % Kernel Results
            FKernel=sumKbeta(KHL , Mu);
            [U6 , ~]= mykernelkmeans(FKernel , CluNum);
            R6 = myNMIACC(U6,Y,CluNum);
            GLMKC_ResK_NLLE(IRate , IBeta , 1) = R6(1);   GLMKC_ResK_NLLE(IRate , IBeta , 2) = R6(2);   GLMKC_ResK_NLLE(IRate , IBeta , 3) = R6(3);
            
            % Ordinary Z Results
            PI = Z > 0;
            Z = Z.*PI;
            [U7] = baseline_spectral_onkernel( abs( (Z + Z') / 2) , CluNum);
            R7 = myNMIACC(U7,Y,CluNum);
            GLMKC_ResZ_NLLE(IRate , IBeta , 1) = R7(1);  GLMKC_ResZ_NLLE(IRate , IBeta , 2) = R7(2);  GLMKC_ResZ_NLLE(IRate , IBeta , 3) = R7(3);
            
        end
    end
    res(:, 3) = [ max(max(GLMKC_ResK_NLLE(: , : , 1))) , max(max(GLMKC_ResK_NLLE(: , : , 2))) ,  max(max(GLMKC_ResK_NLLE(: , : , 3)))];
    res(:, 4) = [ max(max(GLMKC_ResZ_NLLE(: , : , 1))) , max(max(GLMKC_ResZ_NLLE(: , : , 2))) ,  max(max(GLMKC_ResZ_NLLE(: , : , 3)))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%       GLMKC      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  super parameter setting
    [KHL] = LocalKernelCalculationNew(KH , LowRankNum);
    LargestIteration = 15;
    KernelWeight = ones(1 , KerNum) / KerNum;
    avgKer=sumKbeta(KHL,KernelWeight);
    BaseValue = norm(avgKer , 'fro')^2;
    RegularizationValue = 10^(-4) * BaseValue;
    Alpha_T = 2 * BaseValue;
    Beta_T = 2.^(-4 : 4 : 4);
    NNRate = [0.01];
    
    GLMKC_ResZ = zeros(length(NNRate) , length(Beta_T) , 3);
    GLMKC_ResK = zeros(length(NNRate) , length(Beta_T) , 3);
    
    LMKC_ResZ = zeros(length(NNRate) , length(Beta_T) , 3);
    LMKC_ResK = zeros(length(NNRate) , length(Beta_T) , 3);
    
    %% Algorithm
    for IRate = 1 : length(NNRate)
        
        % Preprocessing (Calculating KHL, L , M1)
        CRate = NNRate(IRate);
        [KHL , L , M1] = LGMKCPreprocessFinalNew(KH , LowRankNum);
        
        for IBeta = 1 : length(Beta_T)
            Beta_C = Beta_T(IBeta);
            
            InputKnum = max(2* CluNum , LowRankNum);
            
            % Kernel weight and Z calculation
            [Mu , Z , flag , TotalObj] = LGMKC(KHL, M1, RegularizationValue , Alpha_T , Beta_C , L , LargestIteration , InputKnum);
            
            [Mu2 , Z2 , flag2 , TotalObj2] = LMKC(KHL , M1 ,  RegularizationValue , Alpha_T , Beta_C , L , LargestIteration);
            
            PON = TotalObj(1:end-1) - TotalObj(2 : end);
            DON = sum(PON < -10^(-10)) == 0;
            fprintf('\n\n The Current dataset is : %d; IBeta = %d , IGamma = %d , flag = %d \n\n ' , ICount , IRate , IBeta , DON)
            
            % Result Record
            % Kernel Results
            FKernel=sumKbeta(KHL , Mu);
            [U6 , ~]= mykernelkmeans(FKernel , CluNum);
            R6 = myNMIACC(U6,Y,CluNum);
            GLMKC_ResK(IRate , IBeta , 1) = R6(1);   GLMKC_ResK(IRate , IBeta , 2) = R6(2);   GLMKC_ResK(IRate , IBeta , 3) = R6(3);
            
            % Ordinary Z Results
            PI = Z > 0;
            Z = Z.*PI;
            [U7] = baseline_spectral_onkernel( abs( (Z + Z') / 2) , CluNum);
            R7 = myNMIACC(U7,Y,CluNum);
            GLMKC_ResZ(IRate , IBeta , 1) = R7(1);  GLMKC_ResZ(IRate , IBeta , 2) = R7(2);  GLMKC_ResZ(IRate , IBeta , 3) = R7(3);
            
            % Kernel Results No Krank
            FKernel=sumKbeta(KHL , Mu2);
            [U8 , ~]= mykernelkmeans(FKernel , CluNum);
            R8 = myNMIACC(U8,Y,CluNum);
            LMKC_ResK(IRate , IBeta , 1) = R8(1);   LMKC_ResK(IRate , IBeta , 2) = R8(2);   LMKC_ResK(IRate , IBeta , 3) = R8(3);
            
            % Ordinary Z Results No Krank
            PI = Z2 > 0;
            Z2 = Z2.*PI;
            [U9] = baseline_spectral_onkernel( abs( (Z2 + Z2') / 2) , CluNum);
            R9 = myNMIACC(U9,Y,CluNum);
            LMKC_ResZ(IRate , IBeta , 1) = R9(1);  LMKC_ResZ(IRate , IBeta , 2) = R9(2);  LMKC_ResZ(IRate , IBeta , 3) = R9(3);
        end
    end
    res(:, 5) = [ max(max(LMKC_ResK(: , : , 1))) , max(max(LMKC_ResK(: , : , 2))) ,  max(max(LMKC_ResK(: , : , 3)))];
    res(:, 6) = [ max(max(LMKC_ResZ(: , : , 1))) , max(max(LMKC_ResZ(: , : , 2))) ,  max(max(LMKC_ResZ(: , : , 3)))];
    res(:, 7) = [ max(max(GLMKC_ResK(: , : , 1))) , max(max(GLMKC_ResK(: , : , 2))) ,  max(max(GLMKC_ResK(: , : , 3)))];
    res(:, 8) = [ max(max(GLMKC_ResZ(: , : , 1))) , max(max(GLMKC_ResZ(: , : , 2))) ,  max(max(GLMKC_ResZ(: , : , 3)))];
    
    res = res';
    save(['.\Result\1114_Final_' , dataName , '.mat'], 'res' , 'GLMKC_ResK', 'GLMKC_ResZ', 'LMKC_ResK' ,'LMKC_ResZ', 'Alpha_T', ...
        'Beta_T' , 'GLMKC_ResZ_NLLE' , 'GLMKC_ResK_NLLE' , 'NNRate')
end