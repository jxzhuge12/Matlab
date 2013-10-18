function [d] = aver(a,b)
X = [a+b,13-a,13-b];
n = max(X);
C = rand(n,n);
d = zeros(n,n);
B = zeros(n+2,n+2);
for i = 1:n
    for j = 1:n
        B(i+1,j+1) = C(i,j);
    end
end
for i = 1:n
    for j = 1:n
        d(i,j) = (B(i+1,j+1) + B(i,j+1) + B(i+2,j+1) + B(i+1,j) + B(i+1,j+2))/5;
    end
end