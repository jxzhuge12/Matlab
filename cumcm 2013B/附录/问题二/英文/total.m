function [ansgraph] = total(row)

A = zeros(180,72,209);
for i = 0:2
    for j = 0:9
        for k = 0:9
            if(i * 100 + j * 10 + k <= 208)
                A(:,:,i * 100 + j * 10 + k + 1) = imread([num2str(i),num2str(j),num2str(k)],'bmp');
            end
        end
    end
end
% ¶ÁÈ¡Í¼Ïñ
row = genermatr();
ansgraph = zeros(11,19);
[w,o]=findleft(A);
w = w';
know = 0;
know = input('Do you need to use OCR function?(1/0)');
if(know == 0)
    gr = zeros(180,1);
    tit = input('Where do you want to start?');
    for i = tit:11
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
            for t = 1:209
                row(t,ansgraph(i,d)) = 0;
                row(ansgraph(i,d),t) = 0;
            end
        end
        disp((ansgraph(i,:)));
    end
else
    gr = zeros(180,1);
    tit = input('Where do you want to start?');
    for i = tit:11
        inputnum = 1;
        while(inputnum ~= 0)
            tempj = 1;
            row2 = w(i);
            ansgraph(i,tempj) = row2;
            while(tempj < 19 && ~isempty(find(row(row2,:) ~= 0)))
                p1 = find(row(row2,:) == max(row(row2,:)));
                if(hd (A(:,:,row2),A(:,:,p1)) == 0)
                    row(row2,p1) = 0;
                    continue;
                end
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
            for t = 1:209
                row(t,ansgraph(i,d)) = 0;
                row(ansgraph(i,d),t) = 0;
            end
        end
        disp((ansgraph(i,:)));
    end
end
end