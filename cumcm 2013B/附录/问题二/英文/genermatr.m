function [row] = genermatr()

A = zeros(180,72,209);
row = zeros(209,209);
temprow = 0;
t1=0;

for i = 0:2
    for j = 0:9
        for k = 0:9
            if(i * 100 + j * 10 + k <= 208)
                A(:,:,i * 100 + j * 10 + k + 1) = imread([num2str(i),num2str(j),num2str(k)],'bmp');
            end
        end
    end
end
%读取图像

for i = 1:209
    for j = 1:180
        for k = 1:72
            if(A(j,k,i) < 255)
                A(j,k,i) = 0;
            end
        end
    end
end
%将灰点改成黑点

for i = 1:209
    for j = 1:209
        if(i == j)
            row(i,j) = 0;
        else
            for k = 1:180
                if(A(k,72,i) == 0 && A(k,1,j) == 0)
                    temprow = temprow + 1;
                end
                if(A(k,72,i) == 0) 
                    t1=t1+1; 
                end
                if(A(k,1,j) == 0) 
                    t1=t1+1; 
                end
            end
            if(t1==0) 
                row(i,j)=1;
            else 
                row(i,j) = 2 * temprow / t1;
            end
            t1 = 0;
        end
        temprow = 0;
    end
end
% 对于每两幅图，测量相连的黑点数量

for i = 1:209
    for j = 1:209
        if(i ~= j)
            row(i,j) = row(i,j) + rand()*0.05;
        end
    end
end

[Q]=getblack();
for i = 1:209
    for j = (i+1):209
        if(mod((Q(i) - Q(j)),64) > 2 && 64 - mod((Q(i) - Q(j)),64) > 2)
            row(i,j) = 0;
            row(j,i) = 0;
        end    
    end
end
% 判断两幅图是否可能属于同一行
end