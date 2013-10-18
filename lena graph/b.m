b=zeros(l*l,1);
for i=2:l-1   %方程条件
    for j=2:l-1
b(i+(j-1)*l)=lena(i,j);
  end
end