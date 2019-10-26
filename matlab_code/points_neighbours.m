%求点的邻点。points_nbrs的每一行是点i的所有1领域的点。
function points_nbrs=points_neighbours(face_number,faces,point_number)
degree=points_degree(face_number,faces);
points_degree_temp=zeros(point_number,1);
points_nbrs=zeros(size(point_number,2*max(degree)));
for i=1:face_number
    for j=1:3       
        points_nbrs(faces(i,j),points_degree_temp(faces(i,j))+1)=faces(i,mod(j,3)+1); 
        points_nbrs(faces(i,j),points_degree_temp(faces(i,j))+2)=faces(i,mod(j+1,3)+1);   
        points_degree_temp(faces(i,j))=points_degree_temp(faces(i,j))+2;
    end
end
points_nbrs=sort(points_nbrs,2,'descend');
points_nbrs=points_nbrs(:,1:2:2*max(degree));