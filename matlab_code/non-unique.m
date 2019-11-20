%�ýű�����չʾ�ⲻΨһ�ԡ�ֱ�����м��ɡ�
%�ýű������������һ��ģ�ͣ���ģ��ÿ����ĶȾ�Ϊ4������6���㣬����������׶ƴ�Ӷ�����ͼ�Ρ�
%��Ȼ���ڵ�������ģ������ĵ㲻���棬��ΪҲ������־���ԳƵ����������Ϊ�˱�֤ÿ����Ķȴ��ڵ���4����Ҳ����򵥵�����ˡ�
%�ýű�ʹ�����ϴν����ķ�����⡣���Ǹ÷�����һ������Ч��⣬�ʿ�����Ҫ�໨��ʱ�䡣�����ܵõ�������Ϊ10^-13�Ľ⣬
%������Ϊ�Ѿ��õ��˾��Ƚ⡣������ʾԭʼ��;��Ƚ������ɵ�ģ�͡���ʱ����ģ�ͼ�����ʲ�һ�¡�
%��ˣ���ʹÿ����Ķȴ��ڵ���4�������ĵ㲻���棬�ڳ�����֪������£����Բ�Ψһ��
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
%��������������������������ȡ����ֵ����������������������С�������l2�����Ķ����������С�Ķ��������
[max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(face_number,vertex_dst,vertex_src,points,length);
