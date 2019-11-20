addpath('splsolver')
max_iteration=40;
[points_old,faces,~,~]=readObj('bunny.obj');
length=edge_length(faces,points_old);
%g_is_zero�������е�App D��gֱ����Ϊ0.������40��Ϊ�������Ϊ10^-10����.��ʡ��ʱ����g=0����g_is_zero=true��
g_is_zero=false;
points=BunnyFromMetric(max_iteration,length,faces,g_is_zero);
mesh = Triangular_mesh(faces,length);
%�����ƽ����ת�ԳƳ����������ơ�
mesh.points=fix_rotation(mesh.f_hedge,points_old,mesh.v_dst,mesh.v_src,mesh.v_face,mesh.v_flip,mesh.f_number,points,mesh.faces);
%������
save_obj(mesh.points,mesh.faces,sprintf('BunnyFromMetric_%d.obj',max_iteration));
%��������������������������ȡ����ֵ����������������������С�������l2�����Ķ����������С�Ķ��������
[max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(mesh.f_number,mesh.v_dst,mesh.v_src,mesh.points,mesh.length);


 

