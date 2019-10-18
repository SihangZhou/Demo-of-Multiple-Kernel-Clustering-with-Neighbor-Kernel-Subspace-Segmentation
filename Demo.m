clear
clc
warning off;

DataName = cell(10 , 1);
DataName{1} = 'bbcsport';
DataName{2} = 'bbcsport2view';
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
DataName{13} = 'flower17_DL_fea';

LowRankRate = 0.05;

for ICount = 1
    
    res = zeros(3 , 6);
    
    path = './';
    pathdata = '.';
    addpath(genpath(path));
    dataName = DataName{ICount}; %%% flower17; flower102; CCV; caltech101_mit_numOfClass_10
    %% caltech101_numOfClass_10_Kmatrix
    %% washington texas cornell; wisconsin, UCI_DIGIT
    load([pathdata,'/',dataName,'_Kmatrix.mat'],'KH','Y');
    % load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CluNum = length(unique(Y));
    KerNum = size(KH,3);
    SampleNum = size(KH,1);
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
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%       GLMKC-without local kernel      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %  super parameter setting
%     M1 = MCalculation_New(KH);
%     LargestIteration = 15;
%     KernelWeight = ones(1 , KerNum) / KerNum;
%     avgKer=sumKbeta(KH,KernelWeight);
%     BaseValue = norm(avgKer , 'fro')^2;
%     Alpha_C = 10^(-4) * BaseValue;
%     Beta_T = 0 * BaseValue;
%     Gamma_T = 2.^[-8 , -2 , 2 , 6];
%     NNRate = 0.05;
%     
%     GMKC_ResZ = zeros(length(NNRate) , length(Gamma_T) , 3);
%     GMKC_ResK = zeros(length(NNRate) , length(Gamma_T) , 3);
%     
%     
%     %% Algorithm
%     for IRate = 1 : length(NNRate)
%         for IGamma = 1 : length(Gamma_T)
%             Gamma_C = Gamma_T(IGamma);
%             
%             InputKnum = max(2* CluNum , round(0.05 * SampleNum));
%             
%             % Kernel weight and Z calculation
%             [Mu , Z , flag , TotalObj] = LGMKCNew(KH, M1, Alpha_C , Beta_T , Gamma_C , zeros(size(KH , 1)) , LargestIteration , InputKnum);
%             
%             PON = TotalObj(1:end-1) - TotalObj(2 : end);
%             DON = sum(PON < -10^(-10)) == 0;
%             fprintf('\n\n The Current dataset is : %d; IBeta = %d , IGamma = %d , flag = %d \n\n ' , ICount , IRate , IGamma , DON)
%             
%             % Result Record
%             % Kernel Results
%             FKernel=sumKbeta(KH , Mu);
%             [U14 , ~]= mykernelkmeans(FKernel , CluNum);
%             R14 = myNMIACC(U14,Y,CluNum);
%             GMKC_ResK(IRate , IGamma , 1) = R14(1);   GMKC_ResK(IRate , IGamma , 2) = R14(2);   GMKC_ResK(IRate , IGamma , 3) = R14(3);
%             
%             % Ordinary Z Results
%             PI = Z > 0;
%             Z = Z.*PI;
%             [U15] = baseline_spectral_onkernel( abs( (Z + Z') / 2) , CluNum);
%             R15 = myNMIACC(U15,Y,CluNum);
%             GMKC_ResZ(IRate , IGamma , 1) = R15(1);  GMKC_ResZ(IRate , IGamma , 2) = R15(2);  GMKC_ResZ(IRate , IGamma , 3) = R15(3);
%         end
%     end
%     res(:, 3) = [ max(max(GMKC_ResK(: , : , 1))) , max(max(GMKC_ResK(: , : , 2))) ,  max(max(GMKC_ResK(: , : , 3)))];
%     res(:, 4) = [ max(max(GMKC_ResZ(: , : , 1))) , max(max(GMKC_ResZ(: , : , 2))) ,  max(max(GMKC_ResZ(: , : , 3)))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%       GLMKC     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  super parameter setting
    
    [KHL] = LocalKernelCalculation(KH  , 0.05);
    LargestIteration = 15;
    KernelWeight = ones(1 , KerNum) / KerNum;
    avgKer=sumKbeta(KHL,KernelWeight);
    BaseValue = norm(avgKer , 'fro')^2;
    RegularizationValue = 10^(-4) * BaseValue;
    Alpha_T = 0 * BaseValue;
    Beta_T = 2.^[-8 , -2 , 2 , 6];
    NNRate = [0.01 , 0.03 , 0.09 , 0.11];
    
    %     LowRankRate = RankEstimation(avgKer , RegularizationValue , LowRankRate);
    
    GLMKC_ResZ = zeros(length(NNRate) , length(Beta_T) , 3);
    GLMKC_ResK = zeros(length(NNRate) , length(Beta_T) , 3);
    
    LMKC_ResZ = zeros(length(NNRate) , length(Beta_T) , 3);
    LMKC_ResK = zeros(length(NNRate) , length(Beta_T) , 3);
    
    %% Algorithm
    for IRate = 1 : length(NNRate)
        
        % Preprocessing (Calculating KHL, L , M1)
        CRate = NNRate(IRate);
        [KHL , M1] = LGMKCPreprocessFinal(KH , CRate);
        
        for IBeta = 1 : length(Beta_T)
            Beta_C = Beta_T(IBeta);
            
            InputKnum = max(2* CluNum , round(LowRankRate * SampleNum) );
            
            % Kernel weight and Z calculation
            [Mu , Z , flag , TotalObj] = LGMKCNew(KHL, M1, RegularizationValue , Alpha_T , Beta_C , zeros(SampleNum) , LargestIteration , InputKnum);
            
            [Mu2 , Z2 , flag2 , TotalObj2] = LMKCNew(KHL , M1 ,  RegularizationValue , Alpha_T , Beta_C , zeros(SampleNum) , LargestIteration);
            
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
    %res(:, 7) = [ max(max(GLMKC_ResK(: , : , 1))) , max(max(GLMKC_ResK(: , : , 2))) ,  max(max(GLMKC_ResK(: , : , 3)))];
    %res(:, 8) = [ max(max(GLMKC_ResZ(: , : , 1))) , max(max(GLMKC_ResZ(: , : , 2))) ,  max(max(GLMKC_ResZ(: , : , 3)))];
    
    res = res';
    save(['./Result/19_Final3_' , dataName , '.mat'], 'res' , 'GLMKC_ResK', 'GLMKC_ResZ', 'LMKC_ResK' ,'LMKC_ResZ', 'Alpha_T', ...
        'Beta_T' ,  'NNRate') %'GMKC_ResZ' ,  'GMKC_ResK' ,
end