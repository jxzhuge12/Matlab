function [] = question2chinese()

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
% 读取图像
ansgraph = total(A);
genergraph(A,ansgraph);
finansgraph = total2();
finalgraph = [];
for i = 1:11
    graph = A(:,:,finansgraph(i,1));
    for j = 2:19
        graph = [graph,A(:,:,finansgraph(i,j))];
    end
    finalgraph = [finalgraph;graph];
end

imwrite(finalgraph,'graph.bmp');

end

function [finansgraph] = total2()

A1 = zeros(180,1368,11);
for i = 1:11
    A1(:,:,i) = imread(['graph',num2str(i)],'bmp');
end
col = genermatr2(A1);
finansgraph = zeros(11,1);
gr = zeros(1,1368);

inputnum = 1;
while(inputnum ~= 0)
    tempj = 1;
    col2 = 5;
    finansgraph(tempj) = col2;
    while(tempj < 11 && ~isempty(find(col(col2,:) ~= 0)))
        p1 = find(col(col2,:) == max(col(col2,:)));
        tempj = tempj + 1;
        finansgraph(tempj) = p1;
        col2 = p1;
    end
    graph = A1(:,:,5);
    for l = 2:11
        graph = [graph;gr;A1(:,:,finansgraph(l))];
    end
    imshow(graph);
    inputnum = input('Which picture is wrong? If none, please enter 0.');
    if(inputnum ~= 0)
        col(finansgraph(inputnum),finansgraph(inputnum+1)) = 0;
    end
end

end

function [] = genergraph(A,ansgraph);

for i = 1:11
    graph = A(:,:,ansgraph(i,1));
    for j = 2:19
        graph = [graph,A(:,:,ansgraph(i,j))];
    end
    imwrite(graph,['graph',num2str(i),'.bmp']);
end

end

function [ansgraph] = total(A)

row = genermatr(A);
ansgraph = zeros(11,19);
[w,o]=findleft(A);
w = w';
o = o';
gr = zeros(180,1);
know = input('Where do you want to start?');
for i = know:11
    inputnum = 1;
    while(inputnum ~= 0)
        tempj = 1;
        row2 = w(i);
        ansgraph(i,tempj) = row2;
        while(tempj < 19 && ~isempty(find(row(row2,:) ~= 0)))
            p1 = find(row(row2,:) == max(max(row(row2,:))));
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
            row(ansgraph(i,inputnum),ansgraph(i,inputnum+1)) = 0;
        end
    end
    for d = 1:19
        for t = 1:209
            row(t,ansgraph(i,d)) = 0;
        end
    end
end
end

function [row] = genermatr(A)

row = zeros(209,209);
temprow = 0;
t1=0;
for i = 1:209
    for j = 1:180
        for k = 1:72
            if(A(j,k,i) < 255)
                A(j,k,i) = 0;
            end
        end
    end
end
%将灰点改成黑点
for i = 1:209
    for j = 1:209
        if(i == j)
            row(i,j) = 0;
        else
            for k = 1:180
                if(A(k,72,i) == 0 && A(k,1,j) == 0)
                    temprow = temprow + 1;
                end
                if(A(k,72,i) == 0) 
                    t1=t1+1; 
                end
                if(A(k,1,j) == 0) 
                    t1=t1+1; 
                end
            end
            if(t1==0) 
                row(i,j)=0;
            else 
                row(i,j) = 2 * temprow / t1;
            end
            t1 = 0;
        end
        temprow = 0;
    end
end
for i = 1:209
    for j = 1:209
        row(i,j) = row(i,j) + rand()*0.05;
    end
end
for i = 1:209
    for j = (i+1):209
        [r1,~] = find(A(:,:,i) == 0);
        [r2,~] = find(A(:,:,j) == 0);
        r1 = sort(r1);
        r2 = sort(r2);
        r1 = r1';
        r2 = r2';
        r1 = unique(r1);
        r2 = unique(r2);
        if(length(intersect(r1,r2))/min(length(r1),length(r2)) < 0.9)
            row(i,j) = 0;
            row(j,i) = 0;
        end
    end
end
% 判断两幅图是否可能属于同一行
end

function [col] = genermatr2(A1)

col = zeros(11,11);
tempcol = 0;
temp = 0;
for i = 1:11
    for j = 1:1368
        if (A1(1,j,i) ~= 255) 
            A1(1,j,i) = 0;
        end
        if (A1(180,j,i) ~= 255)
            A1(180,j,i) = 0;
        end
    end
end
for i = 1:11
    for j = 1:11
        if(i == j)
            col(i,j) = 0;
        else
            for h = 1:1368
                if(A1(180,h,i) == 0 && A1(1,h,j) == 0)
                    tempcol = tempcol + 1;
                end
                if(A1(180,h,i) == 0) 
                    temp = temp + 1; 
                end
                if(A1(1,h,j) == 0) 
                    temp = temp + 1; 
                end
            end
            if(temp == 0) 
                col(i,j) = 1;   %两边都无黑点，匹配度为1
            else
                col(i,j) = 2 * tempcol / temp;
            end
        end
        tempcol = 0;
        temp = 0;
    end
end
for i = 1:11
    for j = 1:11
        if(i ~= j)
            col(i,j) = col(i,j) + rand()*0.01;
        end
    end
end
%对于每两幅图，测量相连的黑点数量
end

function [w,o]=findleft(A)
j=1;
x=[];
for i=1:209
    if(A(:,1,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:16
    max=72;
    for j=1:180
        for z=1:72
            if(A(j,z,x(i))==0)
                if(z<max)
                    max=z;
                end
            end
        end
    end
    q(p)=max;
    p=p+1;
end
e=zeros(11,1);
h=zeros(5,1)
e=find(q>10);
h=find(q<10);
w=zeros(11,1);
for i=1:11
    w(i)=x(1,e(1,i));
end
o=zeros(5,1);
for i=1:5
    o(i)=x(1,h(1,i));
end
end
