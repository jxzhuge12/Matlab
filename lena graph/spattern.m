function [rowidx,colidx,entrices] = spattern(A)
%change matrix to sparse matrix
[m,n] = size(A);
rowidx = zeros(1,m+1);
rowidx(1) = 1;
colidx = [];
entrices = [];
number = 1;
for i=1:m
    for j=1:n
        if(A(i,j)~=0 || i == j)
            entrices(number) = A(i,j);
            colidx(number) = j;
            number = number + 1;
        end
    end
    rowidx(i+1) = number;
end
rowidx(m+1)=number-1;