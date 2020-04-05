%该程序是用来测试通过改变对角线长度能否得到满足符合给定边长的图形。
%由于是先随机生成点再计算边长，故满足给定边长的图形总是有的。
%正确率指的是通过dis函数计算得到的error小于10^-5与总次数相比。

%若要进行八面体的实验则用第一行
%faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
faces=[1 2 9;2 3 9;3 4 9;4 5 9;5 6 9;6 7 9;7 8 9;8 1 9;1 2 10;2 3 10;3 4 10;4 5 10;5 6 10;6 7 10;7 8 10;8 1 10];
right=0;
num=1000;%总的次数
%quasi-newton;trust-region
options = optimset('TolFun', 10^-15,'MaxIter',200);
kk1=[];%存储因非实数解而错误的例子
kk2=[];%存储因误差大于10^-5而错误的例子
num1=0;%对应于kk1
num2=0;%对应于kk2
for i=1:num
    %先随机生成点，保证满足对应边长的图形确实存在
    %接下来的两行二选一。第一行是生成随机的十六面体，第二行是生成凸的十六面体。
%     points_sol=rand(10,3);
    z=100;%z用来控制上下两个顶点到x-y平面的距离，距离为rand(1)*z+1
    points_sol=create_convex_Sixteen_face_body(z);%若要对八面体做实验这行可换成points_sol=create_convex_octahedron(z)
    length=edge_length(faces,points_sol);
    %初始化对角线长度，即x0，为全部边长的平均值。
    x0=sum(length(:))/(3*size(length,1));
    f=@(x)dis8(x,points_sol);%若要对八面体做实验这一行可换成f=@(x)dis(x,points_sol);
    
%     x0=norm(points_sol(6,:)-points_sol(5,:));%真实的对角线长度，可选为初值
    [x,fval]=fminunc(f,x0,options);
    %计算输出的error小于10^-5的次数并保证结果是实数
    if fval<10^-5 
        [error,result]=f(x);
        if isreal(result)
            right=right+1;
        else
            kk1=[kk1;points_sol];
            num1=num1+1;
        end
    else
       kk2=[kk2;points_sol]; 
       num2=num2+1;
    end
end
sprintf("正确率为%f%%",right/num*100)

% 
% trimesh(faces, kk1(1:10,1),  kk1(1:10,2),  kk1(1:10,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
% drawnow; pause(0.001);