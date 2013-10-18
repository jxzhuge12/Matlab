function[U]= change(I,n)
U=zeros(n*n,1);

for i=1:n
    U((i-1)*n+1:i*n,1)=I(:,i);
end