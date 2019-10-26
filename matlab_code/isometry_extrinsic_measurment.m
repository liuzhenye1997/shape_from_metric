%求二个相邻平面的二面角。每一行是面上三条半边与相邻面的二面角
function vertex_dihedral=isometry_extrinsic_measurment(vertex_dst,vertex_src,vertex_face,vertex_flip,face_number,points_old,faces)
normal=compute_normal(points_old,faces);

vertex_dihedral=size(3*face_number,1);
vertex_flip=vertex_flip';
vertex_flip=vertex_flip(:);
vertex_face=vertex_face';
vertex_face=vertex_face(:);
for i=1:face_number  
    for j=1:3
        if vertex_flip(3*i+j-3)~=-1
            flip_face=vertex_face(vertex_flip(3*i+j-3));
            N1=normal(vertex_face(3*i+j-3),:);
            N2=normal(flip_face,:);
            P=points_old(vertex_dst(i,j),:)-points_old(vertex_src(i,j),:);
            P_norm=sqrt(sum(P.^2,2));
            E=P./P_norm;
            cross_result=[N2(:,2).*N1(:,3)-N2(:,3).*N1(:,2) N2(:,3).*N1(:,1)-N2(:,1).*N1(:,3) N2(:,1).*N1(:,2)-N2(:,2).*N1(:,1)];
            vertex_dihedral(3*i+j-3)=atan2(sum(E.*cross_result),sum(N1.*N2));
        end
    end
end

