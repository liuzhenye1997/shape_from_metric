%求对偶边的长度。取的是重心，因此将三角形面积三等分，三条边的长度因此为2/3的面积除以对应边的长度
function vertex_dual_length=dual_length_barycentric(length,area)
vertex_dual_length=zeros(size(length));
face_number=size(area,1);
for i=1:face_number
    for j=1:3
        vertex_dual_length(3*i+j-3)=2/3*area(i)/length(i*3+j-3);
    end
end