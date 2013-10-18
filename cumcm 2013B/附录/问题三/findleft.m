function [w]=findleft(A)
j=1;
x=[];
for i=1:418
    if(A(:,1,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:s
    max=72;
    for j=1:180
        for z=1:72
            if(A(j,z,x(i))==0)
                if(z<max)
                    max=z;
                end
            end
        end
    end
    q(p)=max;
    p=p+1;
end
e=zeros(1,22);
[b,n]=sort(q);
e=n(s-21:s);
w=zeros(22,1);
for i=1:22
    w(i)=x(1,e(1,i));
end





    
    
    
        
                
                

