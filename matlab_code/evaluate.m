%��������������������������ȡ����ֵ����������������������С�������l2�����Ķ����������С�Ķ��������
function [max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(face_number,vertex_dst,vertex_src,points,length)
points_last=points;
length_new=sqrt(sum((points_last(vertex_dst',:)-points_last(vertex_src',:)).^2,2));

vertex_log_err=log10(length_new./length);
vertex_rel_err=(length_new./length)-1;
vertex_abs_rel_err=abs(vertex_rel_err);
max_abs_rel_err=max(vertex_abs_rel_err);
min_rel_err=min(vertex_rel_err);
max_rel_err=max(vertex_rel_err);
l2_rel_err=sqrt(sum(vertex_abs_rel_err.^2)/face_number/3);
max_log_err=max(vertex_log_err);
min_log_err=min(vertex_log_err);