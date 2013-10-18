function m-TSPHTW

CityNum=12;		%城市的个数
dislist=[0,667,1600,695,1700,503,1400,952,1500,2219,2200,1500;...
		 667,0,1100,151,1000,2100,2100,2600,1100,1819,1800,952;...
		 1600,1100,0,955,1000,2100,2200,2700,1200,531,550,666;...
		 695,151,955,0,1200,1100,1200,3000,1100,1819,1800,813;...
		 1700,1000,1000,1200,0,2300,2400,2400,1500,1519,1500,781;...
		 503,2100,2100,1100,2300,0,975,547,958,1619,1600,1500;...
		 1400,2100,2200,1200,2400,975,0,454,1000,1719,1700,1600;...
		 952,2600,2700,3000,2400,547,454,0,1500,2119,2100,2000;...
		 1500,1100,1200,1100,1500,958,1000,1500,0,1000,981,297;...
	  	 2219,1819,531,1819,1519,1619,1719,2119,1000,0,19,763;...
		 2200,1800,550,1800,1500,1600,1700,2100,981,19,0,744;...
		 1500,952,666,813,781,1500,1600,2000,297,763,744,0];
			%是距离的12*12的矩阵，距离是从地图上通过确定的路径测出的
Clist=[9.4728 8.5526;...
	   9.2068 8.1005;...
	   9.5651 7.2087;...
	   9.3146 8.03;...
	   8.8817 7.6025;...
       9.8827 8.592;...
	   10.4368 8.3182;...
	   10.3424 8.6656;...
	   9.7185 7.7563;...
	   10.0104 7.2911;...
	   10.0329 7.3291;...
	   9.6982 7.5876];
			%是各个点的坐标

inn=12; 				%初始种群大小
gnmax=500;  			%最大代数
pc=0.8; 				%交叉概率
pm=0.8; 				%变异概率
num=1;				%管理员的个数
minp=2;				%每个管理员管理最少站点的个数
n=CityNum-1;

num_brks=num-1;					% 初始化路线、断点的选择
dof=n-minp*num;          		% 可自由管理的站点数
addto=ones(1,dof+1);
for k=2:num_brks
    addto=cumsum(addto);
end
cum_prob=cumsum(addto)/sum(addto);

sr=zeros(inn,n);        		%初始化种群
sb=zeros(inn,num_brks);  		%初始化断点
for i=1:inn						
    sr(k,:)=randperm(n)+1;		%随机产生单个管理员的路径
    sb(k,:)=randbreaks();		%随机产生断点的位置		
end
[f,p]=objf(s,dislist);

gn=1;
while gn <= gnmax
   for j=1:2:inn
      seln=sel(s,p);  					%选择操作
      scro=cro(s,seln,pc);  				%交叉操作
      scnew(j,:)=scro(1,:);
      scnew(j+1,:)=scro(2,:);
      smnew(j,:)=mut(scnew(j,:),pm);  	%变异操作
      smnew(j+1,:)=mut(scnew(j+1,:),pm);
   end
   s=smnew;  						%产生了新的种群
   [f,p]=objf(s,dislist);  				%计算新种群的适应度
   [fmax,nmax]=max(f);					%记录当前代最好和平均的适应度
   ymean(gn)=1000/mean(f);
   ymax(gn)=1000/fmax;
   x=s(nmax,:);						%记录当前代的最佳个体
   drawTSP(Clist,x,ymax(gn),gn,0); 		%画图
   gn=gn+1;
end
gn=gn-1
end

%------------------------------------------------
%迭代过程中的路径图
function m=drawTSP(Clist,BSF,bsf,p,f)
CityNum=12;
for i=1:CityNum-1
    plot([Clist(BSF(i),1),Clist(BSF(i+1),1)],[Clist(BSF(i),2),Clist(BSF(i+1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g');
plot([Clist(BSF(CityNum),1),Clist(BSF(1),1)],[Clist(BSF(CityNum),2),Clist(BSF(1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g');
title(['12城市TSP']);

hold off;
pause(0.001); %暂停0.001秒，便于观察
end

%------------------------------------------------
%计算当前路径上的损耗
function F=CalDist(dislist,s)

DistanV=0;
n=size(s,2);
for i=1:(n-1)
    DistanV=DistanV+dislist(s(i),s(i+1));
end
ff=1;
tt=0;	%以6点开始计时为例
tmin=[0,0,0,0,0,0,0,0,0,0,0,0];						%到达这些点时间的最小值
tmax=[900,110,900,20,110,900,20,20,110,270,900,900];	%到达这些点时间的最大值
for k=1:num										%分别计算各个管理员的路径消耗
	for i=1:(n-1)
		ff=(tmin(i)<=tt)&(tmax(i)>=tt)&ff;				%判断是否满足给定的时间范围
		tt=tt+3;									%每个点维护时间
		tt=tt+dislist(s(i),s(i+1))/300;					%计算路程上花费的时间
	end
end
DistanV=DistanV+dislist(s(n),s(1));
if ff==0											%不满足时间范围的点
   F=DistanV*100;									%路径消耗上的惩罚
end
ff=ff
F=DistanV;
end

%------------------------------------------------
%计算适应度函数
function [f,p]=objf(s,dislist);

inn=size(s,1);						%读取种群大小。size：获取数组的行数和列数。
for i=1:inn
   f(i)=CalDist(dislist,s(i,:));  		%计算函数值，即适应度
end
f=1000./f';

fsum=0;
for i=1:inn
   fsum=fsum+f(i)^15;				%计算选择概率
end
for i=1:inn
   ps(i)=f(i)^15/fsum;
end

p(1)=ps(1);
for i=2:inn
   p(i)=p(i-1)+ps(i);				%计算累积概率
end
p=p';
end

%--------------------------------------------------
function pcc=pro(pc); 				% pc 交叉概率

test(1:100)=0; 					% 1*100 的零矩阵
l=round(100*pc);					% 100 乘以 pc 然后四舍五入，数值放入l
test(1:l)=1; 						% 第1行第l列的数字改成1
n=round(rand*99)+1; 				% rand产生的是0到1(不包括1)的随机数，所以n的取值是1到100
pcc=test(n); 						% pcc 取值为0或1  
end

%--------------------------------------------------
%“选择”操作
function seln=sel(s,p);

inn=size(p,1);

for i=1:2							%从种群中选择两个个体
   r=rand;  						%产生一个随机数
   prand=p-r;
   j=1;
   while prand(j)<0
       j=j+1;
   end
   seln(i)=j; 						%选中个体的序号
end
end

%------------------------------------------------
%“交叉”操作
function scro=cro(s,seln,pc);

bn=size(s,2);
pcc=pro(pc);  					%根据交叉概率决定是否进行交叉操作，1则是，0则否
scro(1,:)=s(seln(1),:);
scro(2,:)=s(seln(2),:);
if pcc==1
   c1=round(rand*(bn-2))+1;  		%在[1,bn-1]范围内随机产生一个交叉位
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   middle=scro(1,chb1+1:chb2);
   scro(1,chb1+1:chb2)=scro(2,chb1+1:chb2);
   scro(2,chb1+1:chb2)=middle;
   for i=1:chb1
       while find(scro(1,chb1+1:chb2)==scro(1,i))
           zhi=find(scro(1,chb1+1:chb2)==scro(1,i));
           y=scro(2,chb1+zhi);
           scro(1,i)=y;
       end
       while find(scro(2,chb1+1:chb2)==scro(2,i))
           zhi=find(scro(2,chb1+1:chb2)==scro(2,i));
           y=scro(1,chb1+zhi);
           scro(2,i)=y;
       end
   end
   for i=chb2+1:bn
       while find(scro(1,1:chb2)==scro(1,i))
           zhi=find(scro(1,1:chb2)==scro(1,i));
           y=scro(2,zhi);
           scro(1,i)=y;
       end
       while find(scro(2,1:chb2)==scro(2,i))
           zhi=find(scro(2,1:chb2)==scro(2,i));
           y=scro(1,zhi);
           scro(2,i)=y;
       end
   end
end
end

%--------------------------------------------------
%“变异”操作
function snnew=mut(snew,pm);

bn=size(snew,2);
snnew=snew;

pmm=pro(pm);  				%根据变异概率决定是否进行变异操作，1则是，0则否
if pmm==1
   c1=round(rand*(bn-2))+1;  	%在[1,bn-1]范围内随机产生一个变异位
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   x=snew(chb1+1:chb2);
   snnew(chb1+1:chb2)=fliplr(x);
end
end
