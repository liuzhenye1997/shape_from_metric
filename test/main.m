%�ó�������������ͨ���ı�Խ��߳����ܷ�õ�������ϸ����߳���ͼ�Ρ�
%��������������ɵ��ټ���߳�������������߳���ͼ�������еġ�
%��ȷ��ָ����ͨ��dis��������õ���errorС��10^-5���ܴ�����ȡ�
%��errorָ����ͨ����1->2->3��ͨ����1->4->3����·������ĵ�3��λ�õĲ
faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
right=0;
num=10000;%�ܵĴ���
options = optimset('TolFun', 10^-15);
for i=1:num
%     %��������ɵ㣬��֤�����Ӧ�߳���ͼ��ȷʵ����
    points_sol=rand(6,3);
%     points_sol=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 1;0 0 -1];
    length=edge_length(faces,points_sol);
    %��ʼ���Խ��߳��ȣ���x0��Ϊȫ���߳���ƽ��ֵ�ĸ���2����
    x0=sum(length(:))/(3*size(length,1))*sqrt(2);
    f=@(x)dis(x,points_sol);
    [x,fval]=fminunc(f,x0,options);
    %���������errorС��10^-5�Ĵ�������֤�����ʵ��
    if fval<10^-5 && isreal(result)
        [error,result]=f(x);
        right=right+1;
    end
end
sprintf("��ȷ��Ϊ%%%f",right/num*100)
