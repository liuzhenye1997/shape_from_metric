function [resx,resy,resz,points_old,points]=isometry_possion_reconstruction(faces,angle,vertex_prev,vertex_dst,vertex_src,vertex_flip,points,L3,length,face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta)
point_number=size(L3,1);
face_number=size(face_psi1im,1)*size(face_psi1im,2);
vertex_edgevec=construction_edge_direction_per_face(face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta);
vertex_edgevec_flip=vertex_edgevec(vertex_flip,:);
vertex_edgevec=0.5*(vertex_edgevec-vertex_edgevec_flip);
vertex_edgevec=vertex_edgevec.*length;

[points_divx,points_divy,points_divz]=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev);

points_old=points;

LII=L3(2:end,2:end);
points_divxI=points_divx(2:end);
points_divyI=points_divy(2:end);
points_divzI=points_divz(2:end);
solxyzI=LII\[points_divxI points_divyI points_divzI];

solx=[0;solxyzI(:,1)];
soly=[0;solxyzI(:,2)];
solz=[0;solxyzI(:,3)];

resx=L3*solx-points_divx;
resy=L3*soly-points_divy;
resz=L3*solz-points_divz;

points(:,:)=[solx,soly,solz];


area_temp=face_area(points,faces);
weight_temp=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        weight_temp(faces(i,j))=weight_temp(faces(i,j))+area_temp(i)/3;
    end
end

torque_tmp=[0 0 0];
mass_tmp=0;
for i=1:point_number
    torque_tmp=torque_tmp+points(i,:).*weight_temp(i);
    mass_tmp=mass_tmp+weight_temp(i);
end

center_tmp=torque_tmp/mass_tmp;
points=points-center_tmp;

center_tmp=[0 0 0];
torque_tmp=[0 0 0];