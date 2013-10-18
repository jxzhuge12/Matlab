function [gray] = genergray()

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

gray = zeros(209,254);

for i = 1:209
    for j = 1:180
        for k = 1:72
            if(A(j,k,i) ~= 0 && A(j,k,i) ~= 255)
                gray(i,A(j,k,i)) = gray(i,A(j,k,i)) + 1;
            end
        end
    end
end

end

