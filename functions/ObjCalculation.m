function Obj = ObjCalculation(Mu , M1 , M2 , alpha, beta, gamma , Z , L)
Obj = Mu' * (gamma * M1 +  M2) * Mu + alpha * trace(Z' * Z) + beta * trace(Z' * L * Z);
end