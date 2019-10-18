function Kernel = KernelGeneration(K , A)

num = size(K , 1);
Kernel = zeros(num);

for i =1:num
    Ki = zeros(num);
    Ki(A(:,i),A(:,i)) = K(A(:,i),A(:,i));
    Kernel = Kernel + Ki;
end

end