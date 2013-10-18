function [Q]=getblack()
A = zeros(180,72,209);
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
            if(A(j,k,i) ~= 0 && A(j,k,i) ~= 255)
            A(j,k,i) = 0;
            end
        end
    end
end
%将灰点改成黑点
B=zeros(180,1,209);
for i = 1:209
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
C = zeros(180,209);
for i=1:209
    C(:,i) = B(:,1,i);
end

C(76,48) = 20;
C(77,48) = 20;
C(54,65) = 20;

Q=zeros(1,209);
for i = 1:209
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



    
        
            
            
        