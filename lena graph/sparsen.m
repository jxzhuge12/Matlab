function [Amatr] = sparsen(Arowidx,Acolidx,Aentrices)
%稀疏矩阵还原回矩阵
[~,m] = size(Aentrices);
[~,n] = size(Arowidx);
Amatr = zeros(n-1,n-1);
for i = 1:n-1
    if(i == n-1)
        j = Arowidx(i+1) - Arowidx(i) + 1;%A的第i行的个数
    else
        j = Arowidx(i+1) - Arowidx(i);%A的第i行的个数
    end
    for k = 1:j
        Amatr(i,Acolidx(Arowidx(i) - 1 + k)) = Aentrices(Arowidx(i) - 1 + k);
    end
end