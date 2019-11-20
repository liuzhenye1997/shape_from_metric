%This node computes the hedge (houdini vertex) attribute `A` and point attribute `gauss_curvature`. ��
%vertex_A is the argument of the unit complex number r = exp(i*vertex_A) used in the paper .  
%`points_gauss_curvature` is the discrete Gaussian curvature, a.k.a. angle defect, given by 2�� - ��(interior angles) at a point.
function [points_gauss_curvature,vertex_A,vertex_theta]=isometry_levi_civita_connection(angle,faces,point_number,vertex_next,vertex_flip,face_hedge)
%����vertex_theta��������ָ���ǽ������ת��ʵ����������Ƕȣ�Ҳ�ǰ����������нǵ��෴����
vertex_theta=hedge_coordinate(angle,vertex_next,faces,face_hedge);
%����vertex_A
vertex_A=define_parallel_transport_angle_A(vertex_theta,vertex_flip);
vertex_A=fix_pi_pi_ambiguity(vertex_A,vertex_flip);

%������ɢ��˹����
points_gauss_curvature=gauss_curvature(faces,angle,point_number);
points_gauss_curvature=2*pi-points_gauss_curvature;


%% ����vertex��theta��һ���ԡ�������ָ���ǽ������ת��ʵ����������Ƕȣ�Ҳ�ǰ����������нǵ��෴�������е�һ��������������غϡ�������λ��ʵ���Ϸ���
%һ��������������ֵ�ֱ�Ϊ0��A-pi��A+B-2pi��Ȼ��������atan2����ת����-pi��pi֮�䡣����A,B,C�ֱ�Ϊƽ�������ǵĽǶȡ�
function vertex_theta=hedge_coordinate(angle,vertex_next,faces,face_hedge)
face_number=size(faces,1);
vertex_theta=zeros(size(faces))';
he0=face_hedge;
he=he0;
vertex_next_temp=vertex_next';
for i=1:face_number
    cumAngle=0;
    vertex_theta(he(i))=cumAngle;
    he_next=vertex_next_temp(he(i));
    cumAngle=cumAngle+angle(he_next)-pi;
    cumAngle = atan2(sin(cumAngle),cos(cumAngle));
    he(i)=he_next;

    while he(i)~=he0(i)
        vertex_theta(he(i))=cumAngle;
        he_next=vertex_next_temp(he(i));
        cumAngle=cumAngle+angle(he_next)-pi;
        cumAngle = atan2(sin(cumAngle),cos(cumAngle));
        he(i)=he_next;
    end
end
vertex_theta=vertex_theta';

%% ʹvertex_A*vertex_A(vertex_flip)<0�������������ԡ�
function vertex_A=fix_pi_pi_ambiguity(vertex_A,vertex_flip)
vertex_A_temp=vertex_A';
face_number=size(vertex_A,1);
for i=1:face_number
    for j=1:3
        if vertex_flip(3*i+j-3)~=1
            Aopp=vertex_A_temp(vertex_flip(3*(i-1)+j));
            if Aopp*vertex_A(i,j)>0 && 3*i+j-3<vertex_flip(3*i+j-3)
                vertex_A(i,j)=-vertex_A(i,j);
            end
        end
    end
end

%% ��vertex_A,��Ӧ�����е�r(i,j)��
function vertex_A=define_parallel_transport_angle_A(vertex_theta,vertex_flip)
vertex_theta_temp=vertex_theta';
vertex_A=zeros(size(vertex_theta));
face_number=size(vertex_theta,1);
for i=1:face_number
    for j=1:3
        if vertex_flip(3*(i-1)+j)~=-1
            thetaOpp=vertex_theta_temp(vertex_flip(3*(i-1)+j));
            vertex_A(i,j)=thetaOpp-vertex_theta(i,j)-pi;
        end
        vertex_A(i,j)=atan2(sin(vertex_A(i,j)),cos(vertex_A(i,j)));
    end
end

%% ��ÿ����ĸ�˹����
function points_gauss_curvature=gauss_curvature(faces,angle,point_number)
points_gauss_curvature=zeros(1,point_number);
face_number=size(faces,1);
for i=1:face_number
    points_gauss_curvature(faces(i,1))=points_gauss_curvature(faces(i,1))+angle(3*i-2);
    points_gauss_curvature(faces(i,3))=points_gauss_curvature(faces(i,3))+angle(3*i);
    points_gauss_curvature(faces(i,2))=points_gauss_curvature(faces(i,2))+angle(3*i-1);
end