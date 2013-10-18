function [row] = genermatr()

row = zeros(418,418);
temprow = 0;
t1=0;
A = zeros(180,72,418);

for i = 0:2
    for j = 0:9
        for k = 0:9
            if(i * 100 + j * 10 + k <= 208)
                A(:,:,(i * 100 + j * 10 + k + 1) * 2 - 1) = imread([num2str(i),num2str(j),num2str(k),'a'],'bmp');
                A(:,:,(i * 100 + j * 10 + k + 1) * 2) = imread([num2str(i),num2str(j),num2str(k),'b'],'bmp');
            end
        end
    end
end
%读取图像

for i = 1:418
    for j = 1:180
        for k = 1:72
            if(A(j,k,i) ~= 0 && A(j,k,i) ~= 255)
            A(j,k,i) = 0;
            end
        end
    end
end
%将灰点改成黑点

w = findleft(A); % 22个最左
Q = getblack(); % 每幅图的下限

for i = 1:418
    for j = 1:418
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

for i = 1:418
    for j = 1:418
        if(i ~= j)
            row(i,j) = row(i,j) + rand()*0.05;
        end
    end
end

for i = 1:418
    for j = (i+1):418
        if(mod((Q(i) - Q(j)),64) > 2 && 64 - mod((Q(i) - Q(j)),64) > 2)
            row(i,j) = 0;
            row(j,i) = 0;
        end    
    end
end
% 判断两幅图是否可能属于同一行

for i = 1:22
    for j = 1:418
        row(j,w(i)) = 0;
    end
end

end