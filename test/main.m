%�ó�������������ͨ���ı�Խ��߳����ܷ�õ�������ϸ����߳���ͼ�Ρ�
%��������������ɵ��ټ���߳�������������߳���ͼ�������еġ�
%��ȷ��ָ����ͨ��dis��������õ���errorС��10^-5���ܴ�����ȡ�
%��errorָ����ͨ����1->2->3��ͨ����1->4->3����·������ĵ�3��λ�õĲ
%���ڼ���һ���жϽ���Ƿ�Ϊʵ���Ĳ��裬��ȷ�ʴ���½����ӽ�0.
faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
right=0;
num=10000;%�ܵĴ���
options = optimset('TolFun', 10^-15);

kk=0;%�洢���һ�����������
for i=1:num
    %��������ɵ㣬��֤�����Ӧ�߳���ͼ��ȷʵ����
    %�����������ж�ѡһ����һ������������İ����壬�ڶ���������͹�İ����塣
%     points_sol=rand(6,3);
    points_sol=create_convex_octahedron();
    length=edge_length(faces,points_sol);
    %��ʼ���Խ��߳��ȣ���x0��Ϊȫ���߳���ƽ��ֵ��
    x0=sum(length(:))/(3*size(length,1));
%     x0=norm(points_sol(6,:)-points_sol(5,:));%��ʵ�ĶԽ��߳��ȣ���ѡΪ��ֵ
    f=@(x)dis(x,points_sol);
    [x,fval]=fminunc(f,x0,options);
    %���������errorС��10^-5�Ĵ�������֤�����ʵ��
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
sprintf("��ȷ��Ϊ%f%%",right/num*100)

% 
% trimesh(faces, points(:,1), points(:,2), points(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
% drawnow; pause(0.001);