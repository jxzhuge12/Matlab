%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  lena图消噪
%  分别迭代5、50、500、5000次
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = lenagraph()
lena=imread('LenaNoisedSq.jpg');
l = 434;
b=zeros(l*l,1);
for i=2:l-1 %方程条件
    for j=2:l-1
        b(i+(j-1)*l)=lena(i,j);
    end
end
[r,c,e] = matr(434,1,1);
X = gauss(r,c,e,b,5);
X = reshape(X,434,434);
subplot(2,2,1);
imshow(X,[]);
X = gauss(r,c,e,b,50);
X = reshape(X,434,434);
subplot(2,2,2);
imshow(X,[]);
X = gauss(r,c,e,b,500);
X = reshape(X,434,434);
subplot(2,2,3);
imshow(X,[]);
X = gauss(r,c,e,b,5000);
X = reshape(X,434,434);
subplot(2,2,4);
imshow(X,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  五对角矩阵生成（稀疏矩阵三元组格式）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rowidx colidx entrices] = matr(n,t,h)
r=t/(h)^2;
rowidx = zeros(1,n^2+1);
rowidx(1) = 1;
colidx = [];
entrices = [];
temp = 0;
for j = 1:n
    colidx(2*j-1) = j;
    entrices(2*j-1) = r;
    colidx(2*j) = j + n;
    entrices(2*j) = -r;
    rowidx(j+1) = 2 * j + 1;
    temp = temp + 2;
end
for i = 2:n-1
    colidx(temp+1) = (i-1) * n + 1;
    entrices(temp+1) = r;
    colidx(temp+2) = (i-1) * n + 2;
    entrices(temp+2) = -r;
    temp = temp + 2;
    rowidx((i-1)*n+2) = temp + 1;
    for k = 2:n-1
        colidx(temp+1) = (i-2) * n + k;
        entrices(temp+1) = -r;
        colidx(temp+2) = (i-1) * n + k - 1;
        entrices(temp+2) = -r;
        colidx(temp+3) = (i-1) * n + k;
        entrices(temp+3) = 1 + 4 * r;
        colidx(temp+4) = (i-1) * n + k + 1;
        entrices(temp+4) = -r;
        colidx(temp+5) = i * n + k;
        entrices(temp+5) = -r;
        temp = temp + 5;
        rowidx((i-1)*n+k+1) = temp + 1;
    end
    colidx(temp+1) = i * n - 1;
    entrices(temp+1) = -r;
    colidx(temp+2) = i * n;
    entrices(temp+2) = r;
    temp = temp + 2;
    rowidx(i*n+1) = temp + 1;
end
for j = 1:n
    colidx(temp+1) = n^2 - 2 * n + j;
    entrices(temp+1) = r;
    colidx(temp+2) = n^2 - n + j;
    entrices(temp+2) = -r;
    rowidx(n^2-n+j+1) = rowidx(n^2-n+1) + 2 * j;
    temp = temp + 2;
end
rowidx(end) = temp;

%%%%%%%%%%%%
%  高斯迭代
%%%%%%%%%%%%

function [X] = gauss(rowidx,colidx,entrices,b,t)
[~,m] = size(rowidx);
rowidx(m) = rowidx(m) + 1;
X = zeros(1,m-1);
sum = 0;
for j = 1:t %迭代一次
    for i = 1:m-1 %恢复A的第i行
        for k  =rowidx(i):rowidx(i+1)-1
            if(colidx(k) == i)
                a = entrices(k); %用a暂时存储分母A(i,i)的值
                if(a == 0) %避免a=0
                    a = 0.001;
                end
            else
            sum = sum + X(colidx(k)) * entrices(k);
            end
        end
        X(i) = (b(i)-sum)/a;
        sum = 0;
    end
end


