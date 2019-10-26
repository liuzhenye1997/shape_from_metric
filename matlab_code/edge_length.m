%利用点的位置求边长
function length=edge_length(faces,points)
length=zeros(size(faces,1),3);
for i=1:3
    length(:,i)=sqrt(sum((points(faces(:,mod(i-1,3)+1),:)-points(faces(:,mod(i,3)+1),:)).^2,2));
end