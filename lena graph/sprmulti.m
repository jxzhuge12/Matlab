function [Crowidx,Ccolidx,Centrices] = sprmulti(Arowidx,Acolidx,Aentrices,Browidx,Bcolidx,Bentrices)
%n*n稀疏矩阵的乘法A*B（不还原回矩阵）
[~,k] = size(Arowidx);%A和B都是(k-1)阶方阵
[~,l] = size(Bentrices);%B总共有l个元素
Crowidx = zeros(1,k);
Crowidx(1) = 1;
Centrices = [];
Ccolidx = [];
num = 1;
ent = 0;
for i = 1:k-1 %从第1行到第k-1行
    if(i ~= k-1)
        j = Arowidx(i+1) - Arowidx(i);
    else
        j = Arowidx(i+1) - Arowidx(i) + 1;
    end %第i行有j个元素
    for n = 1:k-1 %对B的每一列
        for p = 1:l %B中第p个元素
            if(Bcolidx(p) == n)
                for q = k-1:-1:1
                    if(Browidx(q) <= p)
                        break;
                    end
                end
                for t = 1:j
                    if(q == Acolidx(Arowidx(i)-1+t))
                        ent = ent + Bentrices(p) * Aentrices(Arowidx(i)-1+t);
                    end
                end     
            end
        end
        if(ent ~= 0 || i == n)
            Centrices(num) = ent;
            Ccolidx(num) = n;
            num = num + 1;
        end
        ent = 0;
    end
    Crowidx(i+1) = num;
end
Crowidx(end) = num - 1;