%This node computes the hedge (houdini vertex) attribute `A` and point attribute `gauss_curvature`. 、
%vertex_A is the argument of the unit complex number r = exp(i*vertex_A) used in the paper .  
%`points_gauss_curvature` is the discrete Gaussian curvature, a.k.a. angle defect, given by 2π - Σ(interior angles) at a point.
function [points_gauss_curvature,vertex_A,vertex_theta]=isometry_levi_civita_connection(angle,faces,point_number,vertex_next,vertex_flip,face_hedge)
%计算vertex_theta。该属性指的是将半边旋转到实正半轴所需角度，也是半边与正半轴夹角的相反数。
vertex_theta=hedge_coordinate(angle,vertex_next,faces,face_hedge);
%计算vertex_A
vertex_A=define_parallel_transport_angle_A(vertex_theta,vertex_flip);
vertex_A=fix_pi_pi_ambiguity(vertex_A,vertex_flip);

%计算离散高斯曲率
points_gauss_curvature=gauss_curvature(faces,angle,point_number);
points_gauss_curvature=2*pi-points_gauss_curvature;


%% 计算vertex的theta这一属性。该属性指的是将半边旋转到实正半轴所需角度，也是半边与正半轴夹角的相反数。其中第一条半边与正半轴重合。三角形位于实轴上方。
%一个面的三个顶点的值分别为0，A-pi，A+B-2pi。然后再利用atan2将其转换到-pi到pi之间。其中A,B,C分别为平面三个角的角度。
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

%% 使vertex_A*vertex_A(vertex_flip)<0。消除其歧义性。
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

%% 求vertex_A,对应论文中的r(i,j)。
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

%% 求每个点的高斯曲率
function points_gauss_curvature=gauss_curvature(faces,angle,point_number)
points_gauss_curvature=zeros(1,point_number);
face_number=size(faces,1);
for i=1:face_number
    points_gauss_curvature(faces(i,1))=points_gauss_curvature(faces(i,1))+angle(3*i-2);
    points_gauss_curvature(faces(i,3))=points_gauss_curvature(faces(i,3))+angle(3*i);
    points_gauss_curvature(faces(i,2))=points_gauss_curvature(faces(i,2))+angle(3*i-1);
end