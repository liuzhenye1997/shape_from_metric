function pointvertices=point_to_vertices(max_degree,point_number,faces)
face_number=size(faces,1);
pointvertices=zeros(point_number,max_degree);
degree_temp=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        degree_temp(faces(i,j))=degree_temp(faces(i,j))+1;
        pointvertices(faces(i,j),degree_temp(faces(i,j)))=3*i+j-3;    
    end
end