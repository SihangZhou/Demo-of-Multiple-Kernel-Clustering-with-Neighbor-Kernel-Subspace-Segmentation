function Kernel = NeighborKernelSHNew(CurrentKernel , AVGKer, Num)

NeiNum = Num;

AA = genarateNeighborhood(AVGKer , NeiNum);

Kernel = KernelGeneration(CurrentKernel , AA);

end