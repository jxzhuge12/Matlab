function [Q]=getblack()
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
B=zeros(180,1,418);
for i = 1:418
   for j = 1:180 
        temp=0;
        for k = 1:72
            if(A(j,k,i) == 0)
                temp = temp+1;
            end
        end
        B(j,1,i) = temp;
    end
end
C = zeros(180,418);
for i=1:418
    C(:,i) = B(:,1,i);
end
Q=zeros(1,418);
for i = 1:418
    max=0;
    h=0;
    for j = 1:179
        if(max <(C(j,i)-C(j+1,i)))
            max = C(j,i)-C(j+1,i);
            h=j;
        end
    end
    Q(1,i)=h;
end

    
        
            
            
        