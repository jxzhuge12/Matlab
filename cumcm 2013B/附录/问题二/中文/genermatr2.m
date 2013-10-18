function [col] = genermatr2()

col = zeros(11,11);
tempcol = 0;
temp = 0;

A = zeros(180,1368,11);
for i = 1:11
    A(:,:,i) = imread(['graph',num2str(i)],'bmp');
end

for i = 1:11
    for j = 1:1368
        if (A(1,j,i) ~= 255) 
            A(1,j,i) = 0;
        end
        if (A(180,j,i) ~= 255)
            A(180,j,i) = 0;
        end
    end
end

for i = 1:11
    for j = 1:11
        if(i == j)
            col(i,j) = 0;
        else
            for h = 1:1368
                if(A(180,h,i) == 0 && A(1,h,j) == 0)
                    tempcol = tempcol + 1;
                end
                if(A(180,h,i) == 0) 
                    temp = temp + 1; 
                end
                if(A(1,h,j) == 0) 
                    temp = temp + 1; 
                end
            end
            if(temp == 0) 
                col(i,j) = 0;   %两边都无黑点，匹配度为0
            else
                col(i,j) = 2 * tempcol / temp;
            end
        end
        tempcol = 0;
        temp = 0;
    end
end

for i = 1:11
    for j = 1:11
        if(i ~= j)
            col(i,j) = col(i,j) + rand()*0.01;
        end
    end
end
%对于每两幅图，测量相连的黑点数量

