function [angle,vertex_dual_length,point_area]=isometry_intrinsic_measurement(points,faces,length,vertex_prev,vertex_next,dual_length_type)
angle=vertex_angle(length,vertex_prev,vertex_next);
area=face_area(points,faces);
point_number=size(points,1);
face_number=size(faces,1);
point_area=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        point_area(faces(i,j))=point_area(faces(i,j))+area(i)/3;
    end
end
if strcmp(dual_length_type,'barycentric')
    vertex_dual_length=dual_length_barycentric(length,area);
else
    vertex_dual_length=dual_length_circumcentric(length,area,vertex_next,vertex_prev);
end