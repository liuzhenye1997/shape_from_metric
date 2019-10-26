max_iteration=60;
[points_old,faces,~,~]=readObj('Bunny.obj');
length=edge_length(faces,points_old);
points=BunnyFromMetric(max_iteration,length,faces);
mesh = Triangular_mesh(faces,length);
%将结果平移旋转对称成与输入相似。
mesh.points=fix_rotation(mesh.f_hedge,points_old,mesh.v_dst,mesh.v_src,mesh.v_face,mesh.v_flip,mesh.f_number,points,mesh.faces);
%保存结果
save_obj(mesh.points,mesh.faces,sprintf('BunnyFromMetric_%d.obj',max_iteration));
%测试输出结果，从左到右是相对误差取绝对值后的最大相对误差，最大的相对误差，最小的相对误差，l2误差，最大的对数相对误差，最小的对数相对误差。
[max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(mesh.f_number,mesh.v_dst,mesh.v_src,mesh.points,mesh.length);
%画图
mesh.drawmesh(); 

