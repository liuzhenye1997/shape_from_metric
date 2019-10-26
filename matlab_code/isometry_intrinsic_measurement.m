function [angle,vertex_dual_length,point_area]=isometry_intrinsic_measurement(area,point_number,faces,length,vertex_prev,vertex_next,dual_length_type)
face_number=size(faces,1);
%计算每个面的三个角的角度
angle=vertex_angle(length,vertex_prev,vertex_next);

%计算每个点的面积，定义为其领域三角形之和的1/3。
point_area=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        point_area(faces(i,j))=point_area(faces(i,j))+area(i)/3;
    end
end

%求对偶边的长度。对偶的点选在面的重心或外心
if strcmp(dual_length_type,'barycentric')
    vertex_dual_length=dual_length_barycentric(length,area);
else
    vertex_dual_length=dual_length_circumcentric(length,area,vertex_next,vertex_prev);
end