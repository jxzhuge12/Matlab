function [] = question3()

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
ansgraph = total(row);
genergraph(A,ansgraph);
col = genermatr2();
finansgraph = total2(col);
finalgraph = [];
for i = 1:11
    graph = A(:,:,finansgraph(i,1));
    for j = 2:19
        graph = [graph,A(:,:,finansgraph(i,j))];
    end
    finalgraph = [finalgraph;graph];
end
imwrite(finalgraph,'graph(1).bmp');
for i = 12:22
    graph = A(:,:,finansgraph(i,1));
    for j = 2:19
        graph = [graph,A(:,:,finansgraph(i,j))];
    end
    finalgraph = [finalgraph;graph];
end
imwrite(finalgraph,'graph(2).bmp');
end

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
know = 0;
know = input('Do you need to use OCR function?(1/0)');
if(know == 0)
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
else
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
        for t = 1:418
            row(t,ansgraph(i,d)) = 0;
            row(ansgraph(i,d),t) = 0;
        end
    end
    disp((ansgraph(i,:)));
end
end
end

function [col] = genermatr2()

col = zeros(22,22);
tempcol = 0;
temp = 0;

A = zeros(180,1368,22);
for i = 1:22
    A(:,:,i) = imread(['graph',num2str(i)],'bmp');
end

for i = 1:22
    for j = 1:1368
        if (A(1,j,i) ~= 255) 
            A(1,j,i) = 0;
        end
        if (A(180,j,i) ~= 255)
            A(180,j,i) = 0;
        end
    end
end

for i = 1:22
    for j = 1:22
        if(i == j)
            col(i,j) = 0;
        else
            for h = 1:1368
                if(A(180,h,i) == 0 && A(1,h,j) == 0)
                    tempcol = tempcol + 1;
                end
                if(A(180,h,i) == 0) 
                    temp = temp + 1; 
                end
                if(A(1,h,j) == 0) 
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

for i = 1:22
    for j = 1:22
        if(i ~= j)
            col(i,j) = col(i,j) + rand()*0.01;
        end
    end
end
end

function [row] = genermatr()

row = zeros(418,418);
temprow = 0;
t1=0;
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

w = findleft(A); % 22个最左
Q = getblack(); % 每幅图的下限

for i = 1:418
    for j = 1:418
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
                row(i,j)=1;
            else 
                row(i,j) = 2 * temprow / t1;
            end
            t1 = 0;
        end
        temprow = 0;
    end
end
% 对于每两幅图，测量相连的黑点数量

for i = 1:418
    for j = 1:418
        if(i ~= j)
            row(i,j) = row(i,j) + rand()*0.05;
        end
    end
end

for i = 1:418
    for j = (i+1):418
        if(mod((Q(i) - Q(j)),64) > 2 && 64 - mod((Q(i) - Q(j)),64) > 2)
            row(i,j) = 0;
            row(j,i) = 0;
        end    
    end
end
% 判断两幅图是否可能属于同一行

for i = 1:22
    for j = 1:418
        row(j,w(i)) = 0;
    end
end

end

function [w]=findup(A)
[a,b,c]=size(A);
j=1;
x=[];
for i=1:c
    if(A(1,:,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:s
    max=72;
    for j=1:b
        for z=1:a
            if(A(z,j,x(i))==0)
                if(z<max)
                    max=z;
                end
            end
        end
    end
    q(p)=max;
    p=p+1;
end
e=zeros(1,2);
[d,n]=sort(q);
e=n(s-1:s);
w=zeros(2,1);
for i=1:2
    w(i)=x(1,e(1,i));
end
end

function [w]=findleft(A)
j=1;
x=[];
for i=1:418
    if(A(:,1,i)==255)
    x(j)=i;
    j=j+1;
    end
end
s=length(x);
q=[];
p=1;
for i=1:s
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
e=zeros(1,22);
[b,n]=sort(q);
e=n(s-21:s);
w=zeros(22,1);
for i=1:22
    w(i)=x(1,e(1,i));
end
end

function [Q]=getblack()
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
B=zeros(180,1,418);
for i = 1:418
   for j = 1:180 
        temp=0;
        for k = 1:72
            if(A(j,k,i) == 0)
                temp = temp+1;
            end
        end
        B(j,1,i) = temp;
    end
end
C = zeros(180,418);
for i=1:418
    C(:,i) = B(:,1,i);
end
Q=zeros(1,418);
for i = 1:418
    max=0;
    h=0;
    for j = 1:179
        if(max <(C(j,i)-C(j+1,i)))
            max = C(j,i)-C(j+1,i);
            h=j;
        end
    end
    Q(1,i)=h;
end
end

function [] = genergraph(A,ansgraph)

for i = 1:22
    graph = A(:,:,ansgraph(i,1));
    for j = 2:19
        graph = [graph,A(:,:,ansgraph(i,j))];
    end
    imwrite(graph,['graph',num2str(i),'.bmp']);
end

end

function [startmatr] = gener22matr()

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

w = [28 72 182 212 287 12 19 109 168 273 346 158 179 199 293 332 374 400 8 48 178 229]; % 22个最左
Q = getblack(); % 每幅图的下限

startmatr = zeros(22,22);

for i = 1:22
    for j = 1:22
        if(mod(Q(w(i)) - Q(w(j)),64) < 32)
            startmatr(i,j) = mod(Q(w(i)) - Q(w(j)),64);
            startmatr(j,i) = mod(Q(w(i)) - Q(w(j)),64);
        else
            startmatr(i,j) = 64 - mod(Q(w(i)) - Q(w(j)),64);
            startmatr(i,j) = 64 - mod(Q(w(i)) - Q(w(j)),64);
        end
    end
end

for i = 1:22
    startmatr(i,i) = 32;
end

end

% 输入a1,a2
% 输出两幅图是否能够匹配
function [po] = hd (a1,a2);
 
[h,l1] = size(a1);
[h,l2] = size(a2);
 
aa = [a1,a2];
co = [aa(:,(l1-59):(l1+60))];
 
for i = 1:h
    for j = 1:120
        if(co(i,j) < 64) % 灰度调节
           co(i,j) = 0;
        else co(i,j) = 255;
        end
    end
end
 
edge = zeros(h,1);
ee = 0;
for i = 1:h
    if (co(i,60) + co(i,61) == 0)
        edge(i) = 1;
        ee = ee + 1;
    end
end
 
if (ee > 0)
d = zeros(h,1);
for i = 1:h
    for j = 1:120
        if (co(i,j) == 0)
            d(i) = 1;
            break;
        end
    end
end
 
% 分行操作
 
l = 0;
maxl = 0;
up = 0;
down = 0;
for i =1:h
    if (d(i) == 1) 
        l = l + 1;
    else if (l > maxl)
            maxl = l;
            down = i - 1;
            up = down - l + 1;
        end
         l = 0;
    end
end
ll = floor((up+down)/2) + 32;
if (ll>h) 
    ll=ll-64;
end
 
% 取出边界特征最明显的一行
tl = mod(ll,64);
tt = 0;
tmax = 0;
ttemp = 0;
for i = 1:h
    if (mod(abs(i-tl),64) == 0)
        if (ttemp > tmax)
            tmax = ttemp;
            ttemp = 0;
            tt = i;
        end
    else
        if (edge(i) == 1)
            ttemp = ttemp + 1;
        end
    end
end
 
down = tt;
up = tt - 64;
if (up<=0) 
    up = 1;
end
oo = [co(up:down,:)];    
else
    up = 0;
    down = 0;
    oo = 0;
end
if (oo == 0)
    po=1;
else    
% 取字
left = 60;
right = 61;
[hh,~] = size(oo);
col = zeros(1,120);
 
for i = 1:120
    for j = 1:hh
        if (oo(j,i) == 0)
            col(1,i) = 1;
            break;
        end
    end
end
 
while(col(1,left) == 1)
    left = left - 1;
end
while(col(1,right) == 1)
    right = right + 1;
end
 
ch = [oo(:,left:right)];
[ii,jj] = size(ch);
distance = hausdorff(ch);
    if (distance < 50)
         po = 1;
    else po = 0;
    end
end
end
 
function [currentlow] = hausdorff(i1);
currentlow=1000; 
currentlength=0;
currentnum=1;
[~,length]=size(i1); % 样本长度
for num=1:52    % 对每个模板进行匹配度检测
    imageName=strcat(num2str(num),'.bmp');
    i2=imread(imageName); % 读入模板
    
    [~,testlength]=size(i2); %读入模板长度
    j1=edge(i1,'canny');
    jj1=edge(i2,'canny');
    [aa,cc]=size(j1);
    [aaa,ccc]=size(jj1);
%
for m1=1:aa
    for n1=1:cc
        if j1(m1,n1)>0
            j1(m1,n1)=1;
        else
            j1(m1,n1)=0;
        end   
    end
end
k1=zeros(aa,cc);
k1(1,1)=j1(1,1);
for m2=2:aa-1
    for n2=2:cc-1              
        z=zeros(1,5);
        z(1,1)=j1(m2-1,n2-1)+4;
        z(1,2)=j1(m2-1,n2)+3;
        z(1,3)=j1(m2-1,n2+1)+4;
         z(1,4)=j1(m2,n2-1)+3;
         z(1,5)=j1(m2,n2);
        for x=1:4
           mini=z(1,5); 
            if(mini>z(1,x))
                mini=z(1,x);
            end
        end
         k1(m2,n2)=mini;
    end
end
for mm1=1:aaa
    for nn1=1:ccc
        if jj1(mm1,nn1)>0
            jj1(mm1,nn1)=1;
        else
            jj1(mm1,nn1)=0;
        end   
    end
end
kk1=zeros(aaa,ccc);
kk1(1,1)=jj1(1,1);
for mm2=2:aaa-1
    for nn2=2:ccc-1
        zz=zeros(1,5);
        zz(1,1)=jj1(mm2-1,nn2-1)+4;
        zz(1,2)=jj1(mm2-1,nn2)+3;
        zz(1,3)=jj1(mm2-1,nn2+1)+4;
        zz(1,4)=jj1(mm2,nn2-1)+3;
        zz(1,5)=jj1(mm2,nn2);
        for xx=1:4
           mini2=zz(1,5); 
            if(mini2>zz(1,x))
                mini2=zz(1,x);
            else
                continue
            end   
        end
         kk1(mm2,nn2)=mini2;
    end
end
t=0;
for t1=1:aa
    for t2=1:cc
        if(k1(t1,t2)==1)
            t=t+1;
        end
    end
end
 w=zeros(t,2);
  p=0;
  for t1=1:aa
      for t2=1:cc
            
          if k1(t1,t2)==1
              p=p+1;
              w(p,:)=[t1 t2];
          else
              continue
          end
      end
  end
  tt=0;
for tt1=1:aaa
    for tt2=1:ccc
        if(kk1(tt1,tt2)==1)
            tt=tt+1;
        end
    end
end
 ww=zeros(tt,2);
  pp=0;
  for tt1=1:aaa
      for tt2=1:ccc
            
          if kk1(tt1,tt2)==1
              pp=pp+1;
              ww(pp,:)=[tt1 tt2];
          else
              continue
          end
      end
  end
b=zeros(1,t);
 min_a=((w(1,1)-ww(1,1)).^2+(w(1,2)-ww(1,2)).^2).^(1/2);
for i=1:t
    for j=1:tt
        if ((w(i,1)-ww(j,1)).^2+(w(i,2)-ww(j,2)).^2)<min_a     
            min_a=(w(i,1)-ww(j,1)).^2+(w(i,2)-ww(j,2)).^2;
        end
b(1,i)=min_a;
    end            
end % 求出向量b
max_a=b(1,1); %找出最大值
for i=1:t
    if max_a<b(1,i)
        max_a=b(1,i);    
    end
end
c=zeros(1,tt);
for l=1:tt 
    for k=1:t
        min_b=(ww(1,1)-w(1,1)).^2+(ww(1,2)-w(1,2)).^2;       
        if ((ww(l,1)-w(k,1)).^2+(ww(l,2)-w(k,2)).^2)<min_b
            min_b=(ww(l,1)-w(k,1)).^2+(ww(l,2)-w(k,2)).^2;
        end
c(1,l)=min_b;
    end       
end % 求出向量c
max_b=c(1,1); % 找出最大值
for l=1:tt
    if max_b<c(1,l)
        max_b=c(1,l);
    end
end
if max_a>max_b
    max_hausdorff=max_a;
else
    max_hausdorff=max_b;
end
 
if (max_hausdorff<currentlow)
    currentlow=max_hausdorff;
    currentlength=testlength;
    currentnum=num;
else if (max_hausdorff==currentlow)
        if ((abs(testlength-length))<(abs(currentlength-length)))
            currentlength=testlength;
            currentnum=num;
        end
      end
end
end
end
