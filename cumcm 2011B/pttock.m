function [distance] = pttock(juli,churukou)
    distance = zeros(20,13);
    for i = 1:20 %第i个平台到第j个出口的距离
        for j = 1:13
            [d,path] = floydb(juli,i,churukou(j));
            distance(i,j) = d;
        end
    end
end

