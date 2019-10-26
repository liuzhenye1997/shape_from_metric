%计算每个面的法向量
function normal=compute_normal(points,faces)
edge1=points(faces(:,2),:)-points(faces(:,1),:);
edge2=points(faces(:,3),:)-points(faces(:,1),:);
normal=[edge1(:,2).*edge2(:,3)-edge1(:,3).*edge2(:,2) edge1(:,3).*edge2(:,1)-edge1(:,1).*edge2(:,3) edge1(:,1).*edge2(:,2)-edge1(:,2).*edge2(:,1) ];
normal_length=sqrt(sum(normal.^2,2));
normal=-normal./normal_length;
