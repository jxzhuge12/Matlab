function [sigma] = genersigma(gray)

sigma = zeros(1,209);
miu = zeros(1,209);

for i = 1:209
    sum = 0;
    sum2 = 0;
    sum3 = 0;
    for j = 1:254
        sum = sum + gray(i,j);
        sum2 = sum2 + gray(i,j) * j;
    end
    miu(i) = sum2 / sum;
    for j = 1:254
        sum3 = sum3 + gray(i,j) * (j - miu(i)) ^ 2;
    end
    sigma(i) = miu(i) - sqrt(sum3 / sum);
end

end

