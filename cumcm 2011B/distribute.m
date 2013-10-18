function [res,finnum] = distribute(short,anfa)
    temp = 10000;
    temppos = 1;
    num = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    finnum = [];
    res = zeros(24,20);
    tempres = zeros(24,20);
    times = zeros(1,24);
    fc = 100000;
    sum = zeros(1,24);
    for i = 28:29
        num(21) = i;
        for j = 38:40
            num(22) = j;
            for k = 87:92
                num(24) = k;
                num(23) = 48;
                for l = 1:92
                    for m = 1:24
                        if(short(l,num(m)) < temp)
                            temp = short(l,num(m));
                            temppos = m;
                        end
                    end
                    times(1,temppos) = times(1,temppos) + 1;
                    tempres(temppos,times(1,temppos)) = l;
                    temp = 10000;
                    temppos = 1;
                end %找出了此情况下的分配方案
                for n = 1:24
                    for p = 1:20
                        if(tempres(n,p) ~= 0)
                            sum(1,n) = sum(1,n) + short(num(n),p) * anfa(p);
                        end
                    end
                end
                aver = sum ./ times;
                if(var(aver) < fc)
                    fc = var(aver);
                    res = tempres;
                    finnum = num;
                end
                times = zeros(1,24);
                tempres = zeros(24,20);
                sum = zeros(1,24);
                num(23) = 61;
                for l = 1:92
                    for m = 1:24
                        if(short(l,num(m)) < temp)
                            temp = short(l,num(m));
                            temppos = m;
                        end
                    end
                    times(1,temppos) = times(1,temppos) + 1;
                    tempres(temppos,times(1,temppos)) = l;
                    temp = 10000;
                    temppos = 1;
                end %找出了此情况下的分配方案
                for n = 1:24
                    for p = 1:20
                        if(tempres(n,p) ~= 0)
                            sum(1,n) = sum(1,n) + short(num(n),p) * anfa(p);
                        end
                    end
                end
                aver = sum ./ times;
                if(var(aver) < fc)
                    fc = var(aver);
                    res = tempres;
                    finnum = num;
                end
                times = zeros(1,24);
                tempres = zeros(24,20);
                sum = zeros(1,24);
            end
        end
    end
end

