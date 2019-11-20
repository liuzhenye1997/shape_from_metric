%该脚本是来展示解不唯一性。直接运行即可。
%该脚本先随机建立了一个模型，该模型每个点的度均为4，共有6个点，即两个四棱锥拼接而来的图形。
%当然由于点是随机的，任意四点不共面，因为也不会出现镜像对称的情况。由于为了保证每个点的度大于等于4，这也是最简单的情况了。
%该脚本使用了上次讲述的方法求解。但是该方法不一定能有效求解，故可能需要多花点时间。最终能得到相对误差为10^-13的解，
%可以认为已经得到了精度解。最后会显示原始解和精度解所构成的模型。此时两个模型极大概率不一致。
%因此，即使每个点的度大于等于4，任意四点不共面，在长度已知的情况下，解仍不唯一。
addpath('splsolver')
% points_old=[0 -1 2;0 1 1;1 0 0;-1 0 0;0 -1 0;1 1 -3];
error=1;
iteration=0;
while error>10^-25
    points_old=randn(6,3);
    faces=[1 3 5;1 2 3;1 2 4;1 4 5;2 3 6;2 4 6;3 5 6;4 5 6];
    static=faces(1,:);

    t1=clock();
    point_number=max(faces(:));
    face_number=size(faces,1);

    vertex_src=[faces(:,1);faces(:,2);faces(:,3)];
    vertex_dst=[faces(:,2);faces(:,3);faces(:,1)];

    edge=sort([vertex_src vertex_dst],2);
    [~,index]=sort(edge(:,1)*point_number+edge(:,2));
    edge=edge(index,:);
    edge=edge(1:2:3*face_number,:);
    vertex_src=edge(:,1);
    vertex_dst=edge(:,2);
    length=sqrt(sum((points_old(vertex_src,:)-points_old(vertex_dst,:)).^2,2));
    A=sparse([1:3*face_number/2 1:3*face_number/2],[vertex_dst(:);vertex_src(:)],[ones(3*face_number/2,1);-ones(3*face_number/2,1)],3*face_number/2,point_number);

    A_var=[A(:,1:static(1)-1) A(:,static(1)+1) A(:,static(2)+1:static(3)-1) A(:,static(3)+1:end)];
    A_static=A(:,static);

    solver = splsolver(A_var'*A_var, 'llt');  
    solver.refactorize(A_var'*A_var);

    max_iteration=2*10^4;
    
    points=randn(size(points_old));
    points(static,:)=points_old(static,:);
    points_var=points;
    points_var(static,:)=nan;
    points_static=points_old(static,:);

    while error>10^-25 && iteration<max_iteration
        iteration=iteration+1
        y=A*points;
        y_norm=sqrt(sum(y.^2,2));
        y=y./y_norm.*length;
        y_var=y-A_static*points_static;
        points_var = solver.solve(A_var'*y_var);
        points=[points_var(1:static(1)-1,:);points_static(1,:);points_var(static(1):static(2)-2,:);points_static(2,:);points_var(static(2)-1:static(3)-3,:);points_static(3,:);points_var(static(3)-2:end,:)];
        diff=y-A*points;
        diff=diff(:);
        diff=diff.^2;
        error=sum(diff)
    end
    if error>10^-25
        error=1;
        iteration=0;
    end
end
points_old(:,1)=points_old(:,1)+max(abs(points(:)))+max(abs(points_old(:)));
faces_old=faces+point_number;
faces_sum=[faces;faces_old];
drawmesh(faces_sum, [points;points_old]);
%测试输出结果，从左到右是相对误差取绝对值后的最大相对误差，最大的相对误差，最小的相对误差，l2误差，最大的对数相对误差，最小的对数相对误差。
[max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(face_number,vertex_dst,vertex_src,points,length);
