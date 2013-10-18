% 该附录为问题一的Matlab代码
% 程序运行完，会在该程序的同一目录下生成bmp格式的结果图
function [] = question1()
% 问题一的主函数
graph = zeros(1980,1368);
A = zeros(1980,72,19);
for i = 0:1
    for j = 0:9
        if(i * 10 + j + 1 <= 19)
            A(:,:,i * 10 + j + 1) = imread(['0',num2str(i),num2str(j)],'bmp');
        end
    end
end
% 读取图像
matr = genermatr(A);
[~,ans]=fenpei(matr);
temp = 0;
temprow = 0;
tempcol = 0;
col = 0;
row = 0;
for i = 1:19
    for j = 1:19
        if(matr(i,j) ~= 0)
            temp = 1;
            break;
        end
    end
    if(temp == 0)
        temprow = i;
        break;
    else
        temp = 0;
    end
end
% i是最右的纵切图
for j = 1:19
    for i = 1:19
        if(matr(i,j) ~= 0)
            temp = 1;
            break;
        end
    end
    if(temp == 0)
        tempcol = j;
        break;
    else
        temp = 0;
    end
end
% j是最左的纵切图
col = tempcol;
graph = A(:,:,col);
while(row ~= temprow)
    for i = 1:19
        if(ans(col,i) == 1)
            row = i;
            break;
        end
    end
    graph = [graph,A(:,:,row)];
    col = row;
end
imwrite(graph,'graph.bmp');
end
 
% 生成匹配度矩阵
% 输入19幅纵切图的图形矩阵
% 输出相应的匹配度矩阵
% 匹配度矩阵的第i行第j列表示：第i幅图在左，第j幅图在右的匹配度
function [matr] = genermatr(A)
 
matr = zeros(19,19);
temp = 0;
 
for i = 1:19
    for j = 1:1980
        k = 1;
        if(A(j,k,i) ~= 0 && A(j,k,i) ~= 255)
            A(j,k,i) = 0;
        end
        k = 72;
        if(A(j,k,i) ~= 0 && A(j,k,i) ~= 255)
            A(j,k,i) = 0;
        end
    end
end
% 将灰点改成黑点，增加匹配度
 
for i = 1:19
    for j = 1:19
        if(i == j)
            matr(i,j) = 0;
        else
            for k = 1:1980
                if(A(k,72,i) == 0 && A(k,1,j) == 0)
                    temp = temp + 1;
                end
            end
            matr(i,j) = temp / 1980;
        end
        temp = 0;
    end
end
% 生成匹配度矩阵
 
end
 
% 输入效率矩阵，marix为方阵
% 若效率矩阵中有M，则用一充分大的数代替
% 输出z为最优解，ans为最优分配矩阵
function [z,ans]=fenpei(marix)
 
a=marix;
% 确定矩阵维数
s=length(a);
x=max(max(a));
c=zeros(s,s);
c(:,:)=x;
b=c-a;
a=b;
% 确定矩阵行最小值，进行行减
ml=min(a');
for i=1:s
    a(i,:)=a(i,:)-ml(i);
end
% 确定矩阵列最小值，进行列减
mr=min(a);
for j=1:s
    a(:,j)=a(:,j)-mr(j);
end
 
% start working
num=0;
while(num~=s)  
    index=ones(s);
    index=a&index;
    index=~index;
    flag = zeros(s);
    ans = zeros(s);
    while(sum(sum(index)))
        for i=1:s
            t=0;
            l=0;
            for j=1:s
                if(flag(i,j)==0&&index(i,j)==1)
                    l=l+1;
                    t=j;
                end
            end
            if(l==1)
                flag(:,t)=flag(:,t)+1;
                index(:,t)=0;
                ans(i,t)=1;
            end
        end
        for j=1:s
            t=0;
            r=0;
            for i=1:s
                if(flag(i,j)==0&&index(i,j)==1)
                    r=r+1;
                    t=i;
                end
            end
            if(r==1)
                flag(t,:)=flag(t,:)+1;
                index(t,:)=0;
                ans(t,j)=1;
            end
        end
    end  
    num=sum(sum(ans));
     if(s==num)
        break;
    end
    m=max(max(a));
    for i=1:s
        for j=1:s
            if(flag(i,j)==0)
                if(a(i,j)<m)
                    m=a(i,j);
                end
            end
        end
    end
    for i=1:s
        for j=1:s
            if(flag(i,j)==0)
                a(i,j)=a(i,j)-m;
            end
            if(flag(i,j)==2)
                   a(i,j)=a(i,j)+m;
            end
       end
   end
end
zm=ans.*b;
z=0;
z=sum(sum(zm));
end