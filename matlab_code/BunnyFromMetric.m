function [max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=BunnyFromMetric(max_iteration)

format long
[points,faces]=read_obj('bunny.obj');

face_number=size(faces,1);
point_number=size(points,1);
%create vertex attribute
vertex_face=[1:face_number;1:face_number;1:face_number]';
vertex_src=[faces(:,1),faces(:,3),faces(:,2)];
vertex_dst=[faces(:,3),faces(:,2),faces(:,1)];
vertex_next=[(0:face_number-1)*3+2;(0:face_number-1)*3+3;(0:face_number-1)*3+1]';
vertex_prev=[(0:face_number-1)*3+3;(0:face_number-1)*3+1;(0:face_number-1)*3+2]';
vertex_hedge=[(0:face_number-1)*3+1;(0:face_number-1)*3+2;(0:face_number-1)*3+3]';

vertex_nexthedge=zeros(3*face_number,1);
vertex_dst_temp=vertex_dst';
vertex_src_temp=vertex_src';
vertex_dst_src=[vertex_dst_temp(:),vertex_src_temp(:)];
vertex_dst_src=sort(vertex_dst_src,2);
vertex_value=3*face_number*vertex_dst_src(:,1)+vertex_dst_src(:,2);
[s_vertex_value,index]=sort(vertex_value);
for i=1:3*face_number-1
    if s_vertex_value(i)==s_vertex_value(i+1)
        vertex_nexthedge(index(i))=index(i+1);
        vertex_nexthedge(index(i+1))=index(i);
    end
end

vertex_flip=vertex_nexthedge;


%create face attribute
face_hedge=3*(0:face_number-1)+1;
face_point=vertex_src(face_hedge);
faces_area=face_area(points,faces);

%create point attribute
point_hedge=zeros(point_number);
for i=1:face_number
    for j=1:3
        point_hedge(faces(i,j))=3*(i-1)+j;
    end
end

%calculate hedge length
length=sqrt(sum((points(vertex_dst',:)-points(vertex_src',:)).^2,2));

[~,vertex_dual_length,point_area]=isometry_intrinsic_measurement(points,faces,length,vertex_prev,vertex_next,'circumcentric');

%isometry_build_primal_graph
point_isOld=ones(size(point_number));
point_weight=point_area;

face_dst=vertex_dst;
face_flip=vertex_flip;
face_src=vertex_src;
face_A=0;
face_dihedral=0;
face_theta=0;
face_length=length;
face_dual_length=vertex_dual_length;
face_weight=vertex_dual_length./length;

d0=regular_laplacian(face_number,point_number,vertex_src,vertex_dst);
star0_right=sparse(1:point_number,1:point_number,point_weight);
star1_right=sparse(1:3*face_number,1:3*face_number,face_weight);
L_right=-d0'*(star1_right*d0);
M_right=star0_right;



[angle,vertex_dual_length2,point_area]=isometry_intrinsic_measurement(points,faces,length,vertex_prev,vertex_next,'barycentric');

vertex_theta=hedge_coordinate(angle,vertex_next,faces,face_hedge);

vertex_A=define_parallel_transport_angle_A(vertex_theta,vertex_flip);

vertex_A=fix_pi_pi_ambiguity(vertex_A,vertex_flip);

%计算离散高斯曲率
points_gauss_curvature=gauss_curvature(faces,angle,point_number);
points_gauss_curvature=2*pi-points_gauss_curvature;


points_holonomy=points_gauss_curvature;
points_holonomy=points_holonomy/2;
vertex_A=vertex_A/2;
points_holonomyFromA=zeros(point_number,1);

for i=1:face_number
    for j=1:3
        points_holonomyFromA(vertex_src(i,j))=points_holonomyFromA(vertex_src(i,j))+vertex_A(i,j);
    end
end

points_isBad=label_bad_vertices(points_holonomyFromA,points_holonomy,2);

points_local_isBdy=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        if vertex_flip(3*i+j-3)==-1
            points_local_isBdy(vertex_src(i,j))=1;
            points_local_isBdy(vertex_dst(i,j))=1;
        end
    end
end

root_candidate=find_root_candidate(points_local_isBdy);


%swithch
% root=root_candidate;
% points_distance=-ones(point_number,1);
% points_parent=-ones(point_number,1);
% [points_distance,points_parent]=breadth_first_search(faces,points_parent,points_distance,root);
% max_points_distance=max(points_distance);

coroot=1;
points_co_distance=zeros(point_number,1);
[face_distance,face_parent,face_parenting_vtx,vertex_isCotree]=face_breadth_first_search(vertex_flip,face_number,coroot);

vertex_isCotree=set_isCotree_for_opp_hedges(vertex_isCotree,face_number,vertex_flip);


root=root_candidate;
[points_distance,points_parent,points_parenting_vtx,vertex_isTree]=tree_breadth_first_search(point_number,faces,root,vertex_isCotree);

vertex_generator=zeros(face_number,3);
max_points_distance=max(points_distance);
max_codistance=max(face_distance);

points_cumBad=points_isBad;
for i=max_points_distance:-1:1
    d=i;
    for j=1:point_number
        if points_distance(j)==d
             points_cumBad(points_parent(j))=points_cumBad(points_parent(j))+points_cumBad(j);
        end
    end
end


faces=[faces(:,1),faces(:,3),faces(:,2)];
degree=points_degree(face_number,faces);
max_degree=max(degree);
pointvertices=point_to_vertices(max_degree,point_number,faces);
vertex_flipCount=zeros(face_number*3,1);
vertex_src_temp=vertex_src';
for i=1:point_number
    vtx=1;
    vtx_opp=1;
    for j=1:degree(i)
        vtxtmp=pointvertices(i,j);
        vtxtmp_opp=vertex_flip(vtxtmp);
        if vertex_src_temp(vtxtmp_opp)==points_parent(i)
            vtx=vtxtmp;
            vtx_opp=vtxtmp_opp;
            break;
        end
    end
   
    vertex_flipCount(vtx)=points_cumBad(i);
    vertex_flipCount(vtx_opp)=points_cumBad(i);
end


for i=1:face_number
    for j=1:3
        if mod(vertex_flipCount(3*i+j-3),2)==1
            vertex_A(i,j)=vertex_A(i,j)+pi;
        end        
        vertex_A(i,j)=atan2(sin(vertex_A(i,j)),cos(vertex_A(i,j)));
    end
end

vertex_A_temp=vertex_A';
for i=1:face_number
    for j=1:3
        if vertex_flip(3*i+j-3)~=-1
            Aopp=vertex_A_temp(vertex_flip(3*i+j-3));
            if Aopp*vertex_A(i,j)>0 && 3*i+j-3<vertex_flip(3*i+j-3)
                vertex_A(i,j)=-1*vertex_A(i,j);
            end
        end
    end
end


points_holonomyFromA=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        points_holonomyFromA(vertex_src(i,j))=points_holonomyFromA(vertex_src(i,j))+vertex_A(i,j);
    end
end

evenOrodd=floor((points_holonomy(:)-points_holonomyFromA(:))/pi+0.5);
points_isBad=mod(evenOrodd,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%isometry_spinor1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
points_psi1re=zeros(point_number,1);
points_psi1im=zeros(point_number,1);
points_psi2re=zeros(point_number,1);
points_psi2im=zeros(point_number,1);


points_psi1re(:)=randn(point_number,1);
points_psi1im(:)=randn(point_number,1);
points_psi2re(:)=randn(point_number,1);
points_psi2im(:)=randn(point_number,1);


face_psi1re=zeros(face_number,1);
face_psi1im=zeros(face_number,1);
face_psi2re=zeros(face_number,1);
face_psi2im=zeros(face_number,1);

% 
% face_psi1re(:)=randn(face_number,1);
% face_psi1im(:)=randn(face_number,1);
% face_psi2re(:)=randn(face_number,1);
% face_psi2im(:)=randn(face_number,1);
[face_psi1re,face_psi1im,face_psi2re,face_psi2im]=create_face_psi(0,face_psi1re,face_psi1im,face_psi2re,face_psi2im);
face_psi1re(:)=1:face_number;
face_psi1im(:)=2*(2:face_number++1);
face_psi2re(:)=3*(3:face_number+2);
face_psi2im(:)=4*(4:face_number+3);
face_psi=[face_psi1re face_psi1im face_psi2re face_psi2im];
face_psi_norm=sqrt(sum(face_psi.^2,2));
face_psi1re=face_psi1re./face_psi_norm;
face_psi1im=face_psi1im./face_psi_norm;
face_psi2re=face_psi2re./face_psi_norm;
face_psi2im=face_psi2im./face_psi_norm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%isometry_solver1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[angle,vertex_dual_length,point_area]=isometry_intrinsic_measurement(points,faces,length,vertex_prev,vertex_next,'circumcentric');
[face_dihedral,face_weight,face_dual_length,face_length,face_theta,point_isOld,point_weight,face_dst,face_flip,face_src,face_A]=isometry_build_primal_graph(vertex_dual_length,vertex_src,vertex_flip,vertex_dst,faces_area,point_number,length);
d0=regular_laplacian(face_number,point_number,vertex_src,vertex_dst);
star0=sparse(1:face_number,1:face_number,point_weight);
star1=sparse(1:3*face_number,1:3*face_number,face_weight);
L3=-d0'*(star1*d0);


points_weight=faces_area;
face_src=zeros(face_number,3);
face_dst=zeros(face_number,3);
face_flip=zeros(face_number,3);
face_vtx=zeros(face_number,3);
vertex=zeros(3*face_number,2);
for i=1:face_number
    for j=1:3
        flip=vertex_flip(3*i+j-3);
        dst_prim=(flip-(mod(flip-1,3)+1))/3+1;
        if flip>-1
            face_src(i,j)=i;
            face_dst(i,j)=dst_prim;
            face_vtx(i,j)=3*i+j-3;
            vertex(3*i+j-3,1)=i;
            vertex(3*i+j-3,2)=dst_prim;
        end
    end
end

face_A=vertex_A;
face_dihedral=zeros(face_number,3);
face_theta=vertex_theta;
for i=1:face_number
    for j=1:3
        length_temp=length(3*i+j-3);
        total_dual_length=vertex_dual_length2(3*i+j-3)+vertex_dual_length2(vertex_flip(3*i+j-3));
        face_weight(3*i+j-3)=0.5*length_temp/total_dual_length;
    end
end


face_flip=-ones(face_number,3);
hedges=zeros(face_number,6);
hedges_number=zeros(face_number,1);
for i=1:face_number
    for j=1:3
        hedges_number(face_dst(i,j))=hedges_number(face_dst(i,j))+1;
        hedges(face_dst(i,j),hedges_number(face_dst(i,j)))=3*i+j-3;        
    end
end
for i=1:face_number
    hedges(i,4:6)=3*i-2:3*i;
end
face_dst_temp=face_dst';
for i=1:face_number
    for j=1:3    
        hedge=hedges(face_dst(i,j),:);
        for k=1:6
            if face_dst_temp(hedge(k))==face_src(i,j)
                face_flip(i,j)=hedge(k);             
            end
        end
    end
end


face_ind=1:3*face_number;
[L,M,d0,star1]=connection_laplacian(face_weight,points_weight,face_number,point_number,vertex_src,vertex_dst,vertex_A);
points_psi1re=face_psi1re;
points_psi1im=face_psi1im;
points_psi2re=face_psi2re;
points_psi2im=face_psi2im;

iteration=0;

points_new=zeros(point_number,max_iteration*3);
face_psi=zeros(face_number,max_iteration*4);
while iteration<max_iteration
    if iteration>0
        [face_psi1re,face_psi1im,face_psi2re,face_psi2im,resx,resy,resz,points_old]=solver(iteration,face_hedge,vertex_theta,vertex_face,L3,vertex_flip,vertex_src,vertex_dst,vertex_next,vertex_prev,angle,faces,points_psi1re,points_psi1im,points_psi2re,points_psi2im,face_weight,M,face_ind,star1,face_src,face_dst,face_A,length,points,face_number,face_dihedral,face_theta,face_psi1re,face_psi1im,face_psi2re,face_psi2im);

        points_psi1re=face_psi1re;
        points_psi1im=face_psi1im;
        points_psi2im=face_psi2im;
        points_psi2re=face_psi2re;
        face_psi(:,4*iteration+1:4*iteration+4)=[face_psi1re,face_psi1im,face_psi2re,face_psi2im];
    end

    points_new(:,3*iteration+1:3*iteration+3)=points;
    iteration=iteration+1;
    
    
  
end 

[resx,resy,resz,points_old,points]=isometry_possion_reconstruction(faces,angle,vertex_prev,vertex_dst,vertex_src,vertex_flip,points,L_right,length,face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta);

vertex_olddihedral=isometry_extrinsic_measurment(vertex_dst,vertex_src,vertex_face,vertex_flip,face_number,points_old,faces);
vertex_dihedral=isometry_extrinsic_measurment(vertex_dst,vertex_src,vertex_face,vertex_flip,face_number,points,faces);

vertex_src_temp=vertex_src';
vertex_dst_temp=vertex_dst';
Psrc=points_old(vertex_src_temp(face_hedge),:);
Pdst=points_old(vertex_dst_temp(face_hedge),:);
face_hedge_vec_1=Pdst-Psrc;

vertex_dihedral_prod=vertex_olddihedral.*vertex_dihedral;
dihedral_prod=sum(vertex_dihedral_prod);
reflected=dihedral_prod<0;
points_temp=points;
if reflected==1
    points(:,1)=-points(:,1);
end

normal=compute_normal(points,faces);
vertex_src_temp=vertex_src';
vertex_dst_temp=vertex_dst';
Psrc=points(vertex_src_temp(face_hedge),:);
Pdst=points(vertex_dst_temp(face_hedge),:);
face_hedge_vec_0=Pdst-Psrc;


oldnormal=compute_normal(points_old,faces);
old_face_hedge_vec=face_hedge_vec_1;
Q=compute_dihedral(oldnormal,normal);
QoldE=qrotate(Q,old_face_hedge_vec);
c=sum(batch_normalize(face_hedge_vec_0).*batch_normalize(QoldE),2);
s=sum(batch_normalize(face_hedge_vec_0).*batch_cross(normal,batch_normalize(QoldE)),2);
R=create_quaternion(atan2(s,c),normal);
face_rot=batch_quaternion_multiplication(R,Q);

area_new=face_area(points,faces);

face_rot=face_rot.*area_new;

detail_rot=sum(face_rot,1)/size(face_rot,1);

rot=detail_rot./sqrt(sum(detail_rot.^2));
rot=[rot(1) -rot(2:4)];
points=points_temp;
points=qrotate(rot,points);
if reflected==1
    points(:,1)=-points(:,1);
end


save_obj(points,faces,sprintf('BunnyFromMetric_%d.obj',max_iteration));


points_last=points_new(:,max_iteration*3-2:max_iteration*3);
length_new=sqrt(sum((points_last(vertex_dst',:)-points_last(vertex_src',:)).^2,2));
length=sqrt(sum((points_old(vertex_dst',:)-points_old(vertex_src',:)).^2,2));
vertex_log_err=log10(length_new./length);
vertex_rel_err=(length_new./length)-1;
vertex_abs_rel_err=abs(vertex_rel_err);
max_abs_rel_err=max(vertex_abs_rel_err);
min_rel_err=min(vertex_rel_err);
max_rel_err=max(vertex_rel_err);
l2_rel_err=sqrt(sum(vertex_abs_rel_err.^2)/size(vertex,1));
max_log_err=max(vertex_log_err);
min_log_err=min(vertex_log_err);
