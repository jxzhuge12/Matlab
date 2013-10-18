function [finansgraph] = total2(col)

A1 = zeros(180,1368,22);
for i = 1:22
    A1(:,:,i) = imread(['graph',num2str(i)],'bmp');
end

finansgraph = zeros(22,1);
gr = zeros(1,1368);

inputnum = 1;
while(inputnum ~= 0)
    tempj = 1;
    col2 = 10;
    finansgraph(tempj) = col2;
    while(tempj < 11 && ~isempty(find(col(col2,:) ~= 0)))
        p1 = find(col(col2,:) == max(col(col2,:)));
        tempj = tempj + 1;
        finansgraph(tempj) = p1;
        col2 = p1;
    end
    tempj = 12;
    col2 = 12;
    finansgraph(tempj) = col2;
    while(tempj < 22 && ~isempty(find(col(col2,:) ~= 0)))
        p1 = find(col(col2,:) == max(col(col2,:)));
        tempj = tempj + 1;
        finansgraph(tempj) = p1;
        col2 = p1;
    end
    graph = A1(:,:,10);
    for l = 2:22
        graph = [graph;gr;A1(:,:,finansgraph(l))];
    end
    imshow(graph);
    inputnum = input('Which picture is wrong? If none, please enter 0.');
    if(inputnum ~= 0)
        col(finansgraph(inputnum),finansgraph(inputnum+1)) = 0;
    end
end
end