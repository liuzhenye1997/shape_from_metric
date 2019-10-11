
function vertex_dual_length=dual_length_barycentric(length,area)
vertex_dual_length=zeros(size(length));
face_number=size(area,1);
for i=1:face_number
    for j=1:3
        vertex_dual_length(3*i+j-3)=2/3*area(i)/length(i*3+j-3);
    end
end