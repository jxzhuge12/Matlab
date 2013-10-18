function [c] = nnsprmulti(Arowidx,Acolidx,Aentrices,b)
%A为n*n的矩阵，b为n*1的矩阵。不将A还原回矩阵，直接用稀疏矩阵来与b相乘，得到c
[n,~] = size(b);
c = zeros(n,1);
for i = 1:n
    if(i == n)
        k = Arowidx(i+1) - Arowidx(i) + 1;
    else
        k = Arowidx(i+1) - Arowidx(i);
    end
    for j = 1:k
        c(i,1) = c(i,1) + b(Acolidx(Arowidx(i) - 1 + j)) * Aentrices(Arowidx(i) - 1 + j);
    end
end