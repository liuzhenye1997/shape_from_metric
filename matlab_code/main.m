addpath('splsolver')
max_iteration=40;
name='bunny.obj';
[points_old,faces,~,~]=readObj(name);
length=edge_length(faces,points_old);
%g_is_zero令论文中的App D的g直接设为0.以运行40次为例，误差为10^-10量级.想省点时间设g=0，即g_is_zero=true.  
g_is_zero=false;
if strcmp(name,'torus.obj')
    length=0.05*ones(size(length));
    points=Torus(points_old,max_iteration,length,faces,g_is_zero);
else
    points=BunnyFromMetric(max_iteration,length,faces,g_is_zero);
end
mesh = Triangular_mesh(faces,length);
%将结果平移旋转对称成与输入相似。
mesh.points=fix_rotation(mesh.f_hedge,points_old,mesh.v_dst,mesh.v_src,mesh.v_face,mesh.v_flip,mesh.f_number,points,mesh.faces);
%保存结果
save_obj(mesh.points,mesh.faces,sprintf('BunnyFromMetric_%d.obj',max_iteration));
%测试输出结果，从左到右是相对误差取绝对值后的最大相对误差，最大的相对误差，最小的相对误差，l2误差，最大的对数相对误差，最小的对数相对误差。
[max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(mesh.f_number,mesh.v_dst,mesh.v_src,mesh.points,mesh.length);
fprintf("最大相对误差为%f,L2误差为%f",max_abs_rel_err,l2_rel_err);

 

