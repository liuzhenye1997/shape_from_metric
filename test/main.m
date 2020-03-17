%该程序是用来测试通过改变对角线长度能否得到满足符合给定边长的图形。
%由于是先随机生成点再计算边长，故满足给定边长的图形总是有的。
%正确率指的是通过dis函数计算得到的error小于10^-5与总次数相比。
%而error指的是通过点1->2->3和通过点1->4->3两条路径计算的点3的位置的差。
%由于加了一个判断结果是否为实数的步骤，正确率大幅下降，接近0.
faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
right=0;
num=10000;%总的次数
options = optimset('TolFun', 10^-15);

kk=0;%存储最后一个错误的例子
for i=1:num
    %先随机生成点，保证满足对应边长的图形确实存在
    %接下来的两行二选一。第一行是生成随机的八面体，第二行是生成凸的八面体。
%     points_sol=rand(6,3);
    points_sol=create_convex_octahedron();
    length=edge_length(faces,points_sol);
    %初始化对角线长度，即x0，为全部边长的平均值。
    x0=sum(length(:))/(3*size(length,1));
%     x0=norm(points_sol(6,:)-points_sol(5,:));%真实的对角线长度，可选为初值
    f=@(x)dis(x,points_sol);
    [x,fval]=fminunc(f,x0,options);
    %计算输出的error小于10^-5的次数并保证结果是实数
    if fval<10^-5 
        [error,result]=f(x);
        if isreal(result)
            right=right+1;
        else
            kk=points_sol;
        end
    else
       kk=points_sol; 
    end
end
sprintf("正确率为%f%%",right/num*100)

% 
% trimesh(faces, points(:,1), points(:,2), points(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
% drawnow; pause(0.001);