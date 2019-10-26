max_iteration=60;
[points_old,faces,~,~]=readObj('Bunny.obj');
length=edge_length(faces,points_old);
points=BunnyFromMetric(max_iteration,length,faces);
mesh = Triangular_mesh(faces,length);
%�����ƽ����ת�ԳƳ����������ơ�
mesh.points=fix_rotation(mesh.f_hedge,points_old,mesh.v_dst,mesh.v_src,mesh.v_face,mesh.v_flip,mesh.f_number,points,mesh.faces);
%������
save_obj(mesh.points,mesh.faces,sprintf('BunnyFromMetric_%d.obj',max_iteration));
%��������������������������ȡ����ֵ����������������������С�������l2�����Ķ����������С�Ķ��������
[max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(mesh.f_number,mesh.v_dst,mesh.v_src,mesh.points,mesh.length);
%��ͼ
mesh.drawmesh(); 

