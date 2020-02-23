%该程序是用来测试通过改变对角线长度能否得到满足符合给定边长的图形。
%由于是先随机生成点再计算边长，故满足给定边长的图形总是有的。
%正确率指的是通过dis函数计算得到的error小于10^-5与总次数相比。
%而error指的是通过点1->2->3和通过点1->4->3两条路径计算的点3的位置的差。
faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
right=0;
num=10000;%总的次数
options = optimset('TolFun', 10^-15);
for i=1:num
%     %先随机生成点，保证满足对应边长的图形确实存在
    points_sol=rand(6,3);
%     points_sol=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 1;0 0 -1];
    length=edge_length(faces,points_sol);
    %初始化对角线长度，即x0，为全部边长的平均值的根号2倍。
    x0=sum(length(:))/(3*size(length,1))*sqrt(2);
    f=@(x)dis(x,points_sol);
    [x,fval]=fminunc(f,x0,options);
    %计算输出的error小于10^-5的次数并保证结果是实数
    if fval<10^-5 && isreal(result)
        [error,result]=f(x);
        right=right+1;
    end
end
sprintf("正确率为%%%f",right/num*100)
