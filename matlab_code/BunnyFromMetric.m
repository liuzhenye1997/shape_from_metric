function points=BunnyFromMetric(max_iteration,length,faces,g_is_zero)
%利用边的长度和面创建三角网格
mesh = Triangular_mesh(faces,length);
%求点的邻域面积，对偶边的长度，角度和权重。点的权重为点的邻域面积，面的权重为对偶边长与边长之比。即论文算法2的第一行
[angle,vertex_dual_length,point_area]=isometry_intrinsic_measurement(mesh.f_area,mesh.p_number,mesh.faces,mesh.length,mesh.v_prev,mesh.v_next,'circumcentric');
point_weight=point_area;
face_weight=vertex_dual_length./mesh.length;
%建立levi-civita联络。即论文算法2的第二行
[points_gauss_curvature,vertex_A,vertex_theta]=isometry_levi_civita_connection(angle,mesh.faces,mesh.p_number,mesh.v_next,mesh.v_flip,mesh.f_hedge);
%通过调整τ的正负号使每个vertex star都是figure-0，其中τ对应论文中的τ。即论文算法2的第三四行
[vertex_A,~]=isometry_bundle_sqrt(points_gauss_curvature,vertex_A,mesh.v_src,mesh.v_flip,mesh.p_number,mesh.faces);
%计算laplacian矩阵。
[L_graph,~,~,~]=regular_laplacian(face_weight,point_weight,mesh.f_number,mesh.p_number,mesh.v_src,mesh.v_dst);
%face_src为表面顶点所在面，face_dst为拥有同一条半边相邻的面。face_dihedral为平面三条边与i-j平面的二面角。face_weight为平面的三条半边的权重。
[face_dihedral,~,~,face_dst,face_src]=isometry_build_dual_graph(mesh.length,mesh.f_number,mesh.v_flip,vertex_A,vertex_theta,vertex_dual_length);
face_weight=ones(size(face_weight));
%初始化face_psi，即论文中的λ。即论文算法2的第六行
lambda=isometry_spinor(mesh.f_number);
%迭代使得论文定义的能量变小,输出face_psi,即λ.
iteration=0;
solver1 = 0;
L=0;
index=0;
M=0;
P=0;
while iteration<max_iteration  
    [solver1,lambda,L,index,M,P]=solver(solver1,g_is_zero,mesh,iteration,vertex_theta,L_graph,angle,face_weight,face_src,face_dst,vertex_A,face_dihedral,vertex_theta,lambda,L,index,M,P);
    iteration=iteration+1;
end
%利用求得的face_psi(即λ)
points=isometry_possion_reconstruction(angle,mesh.v_prev,mesh.v_dst,mesh.v_src,mesh.v_flip,L_graph,mesh.length,lambda,mesh.v_face,vertex_theta);






