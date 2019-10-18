function [U] = baseline_spectral_onkernel(K,numClust)

D = diag(sum(K,1));
inv_sqrt_D = sqrt(inv(abs(D)));
L = inv_sqrt_D*K*inv_sqrt_D;
%L = inv(D) * K;
L = (L+L')/2;
%sum(sum(L-L'));
% now do an eigen-decomposition of L
%fprintf('doing eigenvalue decomp...\n');
opts.disp = 0;
[U,V] = eigs(L,numClust,'LA',opts);