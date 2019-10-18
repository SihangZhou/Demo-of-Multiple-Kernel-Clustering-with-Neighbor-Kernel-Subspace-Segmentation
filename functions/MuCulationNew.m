function Mu = MuCulationNew(M1 , M2 , Gamma_C)
% QP problem
opts = optimoptions('quadprog' , 'Algorithm' , 'interior-point-convex' , 'display' , 'off');
% Set up Q
H = M2 + M1 * Gamma_C + eye(size(M2,1))*10^(-9);
% Set up the linear part of the problem
f = [];
A = [];
b = [];
Aeq = ones(1 , size(M2 , 1));
beq = 1;
lb = zeros(size(M2 , 1) , 1);
up = [];
[res] = quadprog(H , f , A , b , Aeq , beq , lb , up , [], opts);
% [Mu] = quadprog(H , f , A , b , Aeq , beq , lb , up);
Mu = res;
end