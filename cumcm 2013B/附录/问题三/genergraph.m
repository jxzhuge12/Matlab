function [] = genergraph(A,ansgraph)

for i = 1:22
    graph = A(:,:,ansgraph(i,1));
    for j = 2:19
        graph = [graph,A(:,:,ansgraph(i,j))];
    end
    imwrite(graph,['graph',num2str(i),'.bmp']);
end

end

