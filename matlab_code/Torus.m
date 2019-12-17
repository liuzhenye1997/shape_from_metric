function points=Torus(points,max_iteration,length,faces,g_is_zero)
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
%����ԭ���ĵ������Ӧ��lambda
lambda=isometry_immersion(length,mesh.f_number,vertex_theta,mesh.v_src,mesh.v_dst,mesh.v_next,mesh.f_hedge,points,mesh.faces);
% ����vertex_A����˵�Ǧ�ʹ����torus��һ��ͬ������
vertex_A=modify_A(mesh,vertex_A,lambda);
%����laplacian����
[L_graph,~,~,~]=regular_laplacian(face_weight,point_weight,mesh.f_number,mesh.p_number,mesh.v_src,mesh.v_dst);
%face_srcΪ���涥�������棬face_dstΪӵ��ͬһ��������ڵ��档face_dihedralΪƽ����������i-jƽ��Ķ���ǡ�face_weightΪƽ���������ߵ�Ȩ�ء�
[face_dihedral,face_weight,~,face_dst,face_src]=isometry_build_dual_graph(mesh.length,mesh.f_number,mesh.v_flip,vertex_A,vertex_theta,vertex_dual_length);
%��ʼ��face_psi���������еĦˡ��������㷨2�ĵ�����
lambda=reshape(1:mesh.f_number*4,mesh.f_number,4);
face_psi_norm=sqrt(sum(lambda.^2,2));
lambda=lambda./face_psi_norm;
%����ʹ�����Ķ����������С,���face_psi,����.
iteration=0;
solver1 = 0;
L=0;
L_index=0;
M=0;
P=0;
AVG=0;
AVG_index=0;
while iteration<max_iteration  
    [solver1,lambda,L,L_index,M,P,AVG,AVG_index]=solver(solver1,g_is_zero,mesh,iteration,vertex_theta,L_graph,angle,face_weight,face_src,face_dst,vertex_A,face_dihedral,vertex_theta,lambda,L,L_index,M,P,AVG,AVG_index);
    iteration=iteration+1;
end
%������õ�face_psi(����)
points=isometry_possion_reconstruction(angle,mesh.v_prev,mesh.v_dst,mesh.v_src,mesh.v_flip,L_graph,mesh.length,lambda,mesh.v_face,vertex_theta);

%% ����vertex_A����˵�Ǧ�ʹ����torus��һ��ͬ������
function vertex_A=modify_A(mesh,vertex_A,lambda)
vertex_A_temp=vertex_A';
B=zeros(size(vertex_A));
nbrs=(mesh.v_flip(:)-1-mod(mesh.v_flip(:)-1,3))/3+1;
for i=1:mesh.f_number
    for j=1:3    
        nbr=nbrs(3*i+j-3);
        IA=[vertex_A_temp(3*i+j-3) 0 0];
        qIA=Quaternion.batch_angleaxis_to_q(-2*IA);
        psi_nbr_transp=Quaternion.multiplication(qIA,lambda(nbr,:));
        qprod=Quaternion.multiplication(Quaternion.qinvert(psi_nbr_transp),lambda(i,:));
        rprod=qprod(1);
        if rprod<0
            vertex_A(i,j)=vertex_A(i,j)+pi;
            B(i,j)=vertex_A(i,j);
            vertex_A(i,j)=atan2(sin(vertex_A(i,j)),cos(vertex_A(i,j)));
        end
    end
end
vertex_A_temp=vertex_A';
for i=1:mesh.f_number
    for j=1:3
        if mesh.v_flip(3*i+j-3)~=-1
            Aopp=vertex_A_temp(mesh.v_flip(3*i+j-3));
            if Aopp*vertex_A(i,j)>0 && 3*i+j-3<mesh.v_flip(3*i+j-3)
                vertex_A(i,j)=-1*vertex_A(i,j);
            end
        end
    end
end



