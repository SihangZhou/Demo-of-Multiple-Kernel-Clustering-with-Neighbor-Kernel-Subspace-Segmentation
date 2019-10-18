function [res]= myNMIACC(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,size(U , 2));
% maxIter = 30;
% tmp1 = zeros(maxIter,1);
% tmp2 = zeros(maxIter,1);
% tmp3 = zeros(maxIter,1);
% for iter = 1:maxIter
%     indx = litekmeans(U_normalized,numclass,'MaxIter',100);
%     indx = indx(:);
%     [newIndx] = bestMap(Y,indx);
%     tmp1(iter) = mean(Y==newIndx);
%     tmp2(iter) = MutualInfo(Y,newIndx);
%     tmp3(iter) = purFuc(Y,newIndx);
% end
% % res = [mean(tmp1), mean(tmp2), mean(tmp3); std(tmp1), std(tmp2), std(tmp3)];
% res = [mean(tmp1), mean(tmp2), mean(tmp3)];
indx = litekmeans(U_normalized,numclass, 'MaxIter', 100, 'Replicates',30);
indx = indx(:);
[newIndx] = bestMap(Y,indx);
res(1) = mean(Y==newIndx);
res(2) = MutualInfo(Y,newIndx);
res(3) = purFuc(Y,newIndx);