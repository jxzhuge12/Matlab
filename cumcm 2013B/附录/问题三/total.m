function [ansgraph] = total(row)

A = zeros(180,72,418);

for i = 0:2
    for j = 0:9
        for k = 0:9
            if(i * 100 + j * 10 + k <= 208)
                A(:,:,(i * 100 + j * 10 + k + 1) * 2 - 1) = imread([num2str(i),num2str(j),num2str(k),'a'],'bmp');
                A(:,:,(i * 100 + j * 10 + k + 1) * 2) = imread([num2str(i),num2str(j),num2str(k),'b'],'bmp');
            end
        end
    end
end
%读取图像

for i = 1:418
    for j = 1:180
        for k = 1:72
            if(A(j,k,i) ~= 0 && A(j,k,i) ~= 255)
            A(j,k,i) = 0;
            end
        end
    end
end
%将灰点改成黑点

ansgraph = zeros(22,19);
w = findleft(A);
w = w';
w(22) = 229;

gr = zeros(180,1);
start = input('Where to start?');
for i = start:22
    inputnum = 1;
    while(inputnum ~= 0)
        tempj = 1;
        row2 = w(i);
        ansgraph(i,tempj) = row2;
        while(tempj < 19 && ~isempty(find(row(row2,:) ~= 0)))
            p1 = find(row(row2,:) == max(row(row2,:)));
            tempj = tempj + 1;
            ansgraph(i,tempj) = p1;
            row2 = p1;
        end
        graph = A(:,:,w(i));
        for l = 2:tempj
            graph = [graph,gr,A(:,:,ansgraph(i,l))];
        end
        imshow(graph);
        inputnum = input('Which picture is wrong? If none, please enter 0.');
        if(inputnum ~= 0)
            if(inputnum == 1000)
                inputnum = input('Which picture is wrong? If none, please enter 0.');
                posit = input('Which picture do you want to show?');
                row(ansgraph(i,inputnum),posit) = 10;
            else
                row(ansgraph(i,inputnum),ansgraph(i,inputnum+1)) = 0;
            end
        end
    end
    for d = 1:19
        for t = 1:418
            row(t,ansgraph(i,d)) = 0;
            row(ansgraph(i,d),t) = 0;
        end
    end
    disp((ansgraph(i,:)));
end

end