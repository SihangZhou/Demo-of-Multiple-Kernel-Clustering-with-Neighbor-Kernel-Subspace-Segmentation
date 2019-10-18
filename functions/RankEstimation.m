function Rank = RankEstimation(Kernel , Para , Rate)

SampleNum = size(Kernel , 1);
KSquare = Kernel*Kernel;
A = KSquare + Para*eye(SampleNum);
AInv = eye(size(A)) / A;
TempN = KSquare * AInv * KSquare;
TempN = (TempN + TempN') / 2;
hehe = eig(TempN);
hehe = sort(real(hehe) , 'descend');
Rank = FindIndex(real(hehe) , Rate);

end

function Index = FindIndex(Vector , Rate)
Leng = length(Vector);
CurrentLow = 1;
CurrentHigh = Leng;

TotalSum = sum(Vector);

flag = 0;
while flag == 0
    CurrentIndex = ceil((CurrentLow + CurrentHigh) / 2) ;
    CurrentRate = sum( Vector( 1 : CurrentIndex ) ) / TotalSum;
    if Rate < CurrentRate
        CurrentHigh = CurrentIndex;
    else
        CurrentLow = CurrentIndex;
    end
    
    flag = abs(CurrentHigh - CurrentLow) <3;
end
Index = CurrentIndex;
end