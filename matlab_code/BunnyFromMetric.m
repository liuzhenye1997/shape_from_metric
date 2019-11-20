function points=BunnyFromMetric(max_iteration,length,faces,g_is_zero)
%���ñߵĳ��Ⱥ��洴����������
mesh = Triangular_mesh(faces,length);
%���������������ż�ߵĳ��ȣ��ǶȺ�Ȩ�ء����Ȩ��Ϊ���������������Ȩ��Ϊ��ż�߳���߳�֮�ȡ��������㷨2�ĵ�һ��
[angle,vertex_dual_length,point_area]=isometry_intrinsic_measurement(mesh.f_area,mesh.p_number,mesh.faces,mesh.length,mesh.v_prev,mesh.v_next,'circumcentric');
point_weight=point_area;
face_weight=vertex_dual_length./mesh.length;
%����levi-civita���硣�������㷨2�ĵڶ���
[points_gauss_curvature,vertex_A,vertex_theta]=isometry_levi_civita_connection(angle,mesh.faces,mesh.p_number,mesh.v_next,mesh.v_flip,mesh.f_hedge);
%ͨ�������ӵ�������ʹÿ��vertex star����figure-0�����ЦӶ�Ӧ�����еĦӡ��������㷨2�ĵ�������
[vertex_A,~]=isometry_bundle_sqrt(points_gauss_curvature,vertex_A,mesh.v_src,mesh.v_flip,mesh.p_number,mesh.faces);
%����laplacian����
[L_graph,~,~,~]=regular_laplacian(face_weight,point_weight,mesh.f_number,mesh.p_number,mesh.v_src,mesh.v_dst);
%face_srcΪ���涥�������棬face_dstΪӵ��ͬһ��������ڵ��档face_dihedralΪƽ����������i-jƽ��Ķ���ǡ�face_weightΪƽ���������ߵ�Ȩ�ء�
[face_dihedral,~,~,face_dst,face_src]=isometry_build_dual_graph(mesh.length,mesh.f_number,mesh.v_flip,vertex_A,vertex_theta,vertex_dual_length);
face_weight=ones(size(face_weight));
%��ʼ��face_psi���������еĦˡ��������㷨2�ĵ�����
lambda=isometry_spinor(mesh.f_number);
%����ʹ�����Ķ����������С,���face_psi,����.
iteration=0;
while iteration<max_iteration  
    lambda=solver(g_is_zero,mesh,iteration,vertex_theta,L_graph,angle,face_weight,face_src,face_dst,vertex_A,face_dihedral,vertex_theta,lambda);
    iteration=iteration+1;
end
%������õ�face_psi(����)
points=isometry_possion_reconstruction(angle,mesh.v_prev,mesh.v_dst,mesh.v_src,mesh.v_flip,L_graph,mesh.length,lambda,mesh.v_face,vertex_theta);






