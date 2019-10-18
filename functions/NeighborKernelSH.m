function Kernel = NeighborKernelSH(CurrentKernel , AVGKer, Rate)

NeiNum = round( Rate * size(CurrentKernel , 1) );

AA = genarateNeighborhood(AVGKer , NeiNum);

Kernel = KernelGeneration(CurrentKernel , AA);

end