%�ó�������������ͨ���ı�Խ��߳����ܷ�õ�������ϸ����߳���ͼ�Ρ�
%��������������ɵ��ټ���߳�������������߳���ͼ�������еġ�
%���������1000��Ҳ�޷��õ�һ��������ϸ����߳���ͼ�Ρ�����ȷ��Ϊ0.
faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
right=0;
num=1000;
for i=1:num
    points_sol=rand(6,3);
%     points_sol=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 1;0 0 -1];
    length=edge_length(faces,points_sol);
    x0=sum(length(:))/(3*size(length,1));
    f=@(x)dis(x,points_sol);
    [x,fval]=fminunc(f,x0);
    if fval<10^-5
        [error,result]=f(x);
        right=right+1;
    else
    end
end
sprintf("��ȷ��Ϊ%%%f",right/num*100)
