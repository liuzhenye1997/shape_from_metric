%将输出结果旋转对称到与输入相似的角度和位置
function points=fix_rotation(face_hedge,points_old,vertex_dst,vertex_src,vertex_face,vertex_flip,face_number,points,faces)
point_number=size(points,1);
length_new=zeros(size(faces,1),3);
for i=1:3
    length_new(:,i)=sqrt(sum((points(faces(:,mod(i,3)+1),:)-points(faces(:,mod(i+1,3)+1),:)).^2,2));
end
area=face_area(length_new);
point_weight=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        point_weight(faces(i,j))=point_weight(faces(i,j))+area(i)/3;
    end
end

torque_tmp=[0 0 0];
mass_tmp=0;
for i=1:point_number
    torque_tmp=torque_tmp+points(i,:).*point_weight(i);
    mass_tmp=mass_tmp+point_weight(i);
end

center_tmp=torque_tmp/mass_tmp;
points=points-center_tmp;

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
Q=Quaternion.compute_dihedral(oldnormal,normal);
QoldE=Quaternion.qrotate(Q,old_face_hedge_vec);
c=sum(batch_normalize(face_hedge_vec_0).*batch_normalize(QoldE),2);
s=sum(batch_normalize(face_hedge_vec_0).*batch_cross(normal,batch_normalize(QoldE)),2);
R=Quaternion.create_quaternion(atan2(s,c),normal);
face_rot=Quaternion.batch_quaternion_multiplication(R,Q);

length=edge_length(faces,points);
area_new=face_area(length);

face_rot=face_rot.*area_new;

detail_rot=sum(face_rot,1)/size(face_rot,1);

rot=detail_rot./sqrt(sum(detail_rot.^2));
rot=[rot(1) -rot(2:4)];
points=points_temp;
points=Quaternion.qrotate(rot,points);
if reflected==1
    points(:,1)=-points(:,1);
end