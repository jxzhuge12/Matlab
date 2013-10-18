function [] = draw(ans)

A = zeros(180,1368,11);
for i = 1:11
    A(:,:,i) = imread(['graph',num2str(i)],'bmp');
end

row = 2;
graph = A(:,:,row);
for i = 2:11
    for j = 1:11
        if(ans(row,j) == 1)
            row = j;
            break;
        end
    end
    graph = [graph;A(:,:,row)];
end
imshow(graph);

end

