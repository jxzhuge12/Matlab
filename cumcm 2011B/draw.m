function [] = draw(quanzuobiao,luxian)
    [m,~] = size(luxian);
    [n,~] = size(quanzuobiao);
    X = quanzuobiao(:,1);
    Y = quanzuobiao(:,2);
    for i = 1:n
        plot(X,Y,'x');hold on;
        text(X(i),Y(i),num2str(i));
    end
    for i = 1:m
        plot([quanzuobiao(luxian(i,1),1),quanzuobiao(luxian(i,2),1)],[quanzuobiao(luxian(i,1),2),quanzuobiao(luxian(i,2),2)]);hold on;
    end
end

