function [w c e]= matrixsb(n,t,h)
A=sparse(zeros(n,n));
r=t/(h)^2;
for i=1:n
    A(i,i)=r;
    A(i,n+i)=-r;
end
for j=1:n-2
    A(n*j+1,n*j+1)=r;
    A(n*j+1,n*j+2)=-r;
    for i=n*j+2:n*j+4
        A(i,i-5)=-r;
        A(i,i-1)=-r;
        A(i,i)=1+4*r;
        A(i,i+1)=-r;
        A(i,i+5)=-r;
    end
    A(n*j+5,n*j+4)=-r;
    A(n*j+5,n*j+5)=r;
end
for i=n*(n-1):n*n
    A(i,i)=-r;
    A(i,i-n)=r;
end
[w c e]=spattern(A);