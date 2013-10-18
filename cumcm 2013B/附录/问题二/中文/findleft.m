function [w,o]=findleft(A)
j=1;
x=[];
for i=1:209
    if(A(:,1,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:16
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
e=zeros(11,1);
h=zeros(5,1)
e=find(q>10);
h=find(q<10);
w=zeros(11,1);
for i=1:11
    w(i)=x(1,e(1,i));
end
o=zeros(5,1);
for i=1:5
    o(i)=x(1,h(1,i));
end




    
    
    
        
                
                

