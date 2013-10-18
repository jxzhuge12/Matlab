function [w]=findup(A)
[a,b,c]=size(A);
j=1;
x=[];
for i=1:c
    if(A(1,:,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:s
    max=72;
    for j=1:b
        for z=1:a
            if(A(z,j,x(i))==0)
                if(z<max)
                    max=z;
                end
            end
        end
    end
    q(p)=max;
    p=p+1;
end
e=zeros(1,2);
[d,n]=sort(q);
e=n(s-1:s);
w=zeros(2,1);
for i=1:2
    w(i)=x(1,e(1,i));
end