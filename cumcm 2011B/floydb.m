function [d,path]=floydb(blinjie,sp,ep) %path 是路由矩阵 %path 是路由矩阵
    n=size(blinjie,1);
    D=blinjie;
    path=zeros(n,n);
    for i=1:n
        for j=1:n
            if D(i,j)~=inf
                path(i,j)=j;
            end
        end
    end
    for k=1:n
        for i=1:n
            for j=1:n
                if D(i,k)+D(k,j)<D(i,j)
                    D(i,j)=D(i,k)+D(k,j);
                    path(i,j)=path(i,k);
                end
            end
        end
    end
    p =[sp];
    mp =sp;
    for k=1:n
        if mp ~= ep
            d=path(mp,ep);
            p=[p,d];
            mp=d;
        end
    end
    d = D(sp,ep);
    path = p;
end