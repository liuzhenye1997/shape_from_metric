%�ó�������������ͨ���ı�Խ��߳����ܷ�õ�������ϸ����߳���ͼ�Ρ�
%��������������ɵ��ټ���߳�������������߳���ͼ�������еġ�
%��ȷ��ָ����ͨ��dis��������õ���errorС��10^-5���ܴ�����ȡ�

%��Ҫ���а������ʵ�����õ�һ��
%faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
faces=[1 2 9;2 3 9;3 4 9;4 5 9;5 6 9;6 7 9;7 8 9;8 1 9;1 2 10;2 3 10;3 4 10;4 5 10;5 6 10;6 7 10;7 8 10;8 1 10];
right=0;
num=1000;%�ܵĴ���
%quasi-newton;trust-region
options = optimset('TolFun', 10^-15,'MaxIter',200);
kk1=[];%�洢���ʵ��������������
kk2=[];%�洢��������10^-5�����������
num1=0;%��Ӧ��kk1
num2=0;%��Ӧ��kk2
for i=1:num
    %��������ɵ㣬��֤�����Ӧ�߳���ͼ��ȷʵ����
    %�����������ж�ѡһ����һ�������������ʮ�����壬�ڶ���������͹��ʮ�����塣
%     points_sol=rand(10,3);
    z=100;%z�������������������㵽x-yƽ��ľ��룬����Ϊrand(1)*z+1
    points_sol=create_convex_Sixteen_face_body(z);%��Ҫ�԰�������ʵ�����пɻ���points_sol=create_convex_octahedron(z)
    length=edge_length(faces,points_sol);
    %��ʼ���Խ��߳��ȣ���x0��Ϊȫ���߳���ƽ��ֵ��
    x0=sum(length(:))/(3*size(length,1));
    f=@(x)dis8(x,points_sol);%��Ҫ�԰�������ʵ����һ�пɻ���f=@(x)dis(x,points_sol);
    
%     x0=norm(points_sol(6,:)-points_sol(5,:));%��ʵ�ĶԽ��߳��ȣ���ѡΪ��ֵ
    [x,fval]=fminunc(f,x0,options);
    %���������errorС��10^-5�Ĵ�������֤�����ʵ��
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
sprintf("��ȷ��Ϊ%f%%",right/num*100)

% 
% trimesh(faces, kk1(1:10,1),  kk1(1:10,2),  kk1(1:10,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
% drawnow; pause(0.001);