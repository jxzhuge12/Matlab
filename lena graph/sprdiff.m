function [Crowidx,Ccolidx,Centrices] = sprdiff(Arowidx,Acolidx,Aentrices,Browidx,Bcolidx,Bentrices)
% C = A - B
[~,k] = size(Arowidx);
[~,r] = size(Aentrices);
[~,s] = size(Bentrices);
Crowidx = zeros(1,k);
Crowidx(1) = 1;
Ccolidx = [];
Centrices = [];
num = 1;
for i = 1:k-1
    if(i ~= k-1)
        a = Arowidx(i+1) - Arowidx(i);%A的第i行个数
        b = Browidx(i+1) - Browidx(i);%B的第i行个数
    else
        a = Arowidx(i+1) - Arowidx(i)+1;%A的第i行个数
        b = Browidx(i+1) - Browidx(i)+1;%B的第i行个数
    end
    Aj = 1;
    Bj = 1;
    while(Aj <= a && Bj <= b)
        if(Acolidx(Arowidx(i)-1+Aj) < Bcolidx(Browidx(i)-1+Bj))
            Ccolidx(num) = Acolidx(Arowidx(i)-1+Aj);
            Centrices(num) = Aentrices(Arowidx(i)-1+Aj);
            num = num + 1;
            Aj = Aj + 1;
        elseif(Acolidx(Arowidx(i)-1+Aj) > Bcolidx(Browidx(i)-1+Bj))
            Ccolidx(num) = Bcolidx(Browidx(i)-1+Bj);
            Centrices(num) = -Bentrices(Browidx(i)-1+Bj);
            num = num + 1;
            Bj = Bj + 1;
        else
            if((Aentrices(Arowidx(i)-1+Aj) - Bentrices(Browidx(i)-1+Bj) ~= 0) || (Acolidx(Arowidx(i)-1+Aj) == i))
                Ccolidx(num) = Acolidx(Arowidx(i)-1+Aj);
                Centrices(num) = Aentrices(Arowidx(i)-1+Aj) - Bentrices(Browidx(i)-1+Bj);
                num = num + 1;
                Aj = Aj + 1;
                Bj = Bj + 1;
            else
                Aj = Aj + 1;
                Bj = Bj + 1;
            end
        end
    end
    if(Aj > a)
        for z = Bj:b
            Ccolidx(num) = Bcolidx(Browidx(i)-1+z);
            Centrices(num) = -Bentrices(Browidx(i)-1+z);
            num = num + 1;
        end
    else
        for z = Aj:a
            Ccolidx(num) = Acolidx(Arowidx(i)-1+z);
            Centrices(num) = Aentrices(Arowidx(i)-1+z);
            num = num + 1;
        end
    end
    Crowidx(i+1) = num;
end
Crowidx(k) = num - 1;
