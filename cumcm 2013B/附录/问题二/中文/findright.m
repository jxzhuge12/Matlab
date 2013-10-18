function [t,l]=findright(A)
j=1;
x=[];
for i=1:209
    if(A(:,72,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:15
    max=1;
    for j=1:180
        for z=1:72
            if(A(j,z,x(i))==0)
                if(z>max)
                    max=z;
                end
            end
        end
    end
    q(p)=max;
    p=p+1;
end
e=zeros(11,1);
e=find(q<66);
h=zeros(4,1);
h=find(q>66);
t=zeros(1,11);
for i=1:11
    t(i)=x(1,e(1,i));
end
l=zeros(1,4);
for i=1:4
    l(i)=x(1,h(1,i));
end