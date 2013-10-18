function [res] = res(juli)
    temp = 10000;
    temppos = 1;
    res = zeros(20,92);
    times = zeros(1,20);
    for i = 1:92
        for j = 1:20
            [d,~]=floydb(juli,i,j);
            if(d < temp)
                temp = d;
                temppos = j;
            end
        end
        times(1,temppos) = times(1,temppos) + 1;
        res(temppos,times(1,temppos)) = i;
        temp = 10000;
        temppos = 1;
    end
end

