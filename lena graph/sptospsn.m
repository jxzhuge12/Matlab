function [rowidx colidx entrices] = sptospsn(row,col,ent)
[~,k] = size(row);
n = row(k);
rowidx = zeros(1,n+1);
rowidx(1) = 1;
tem = 1;
colidx = col;
entrices = ent;
for i = 2:k
    if(row(i-1) ~= row(i))
        tem = tem + 1;
        rowidx(tem) = i;
    end
end
rowidx(n+1) = k;
