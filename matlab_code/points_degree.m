function degree=points_degree(face_number,faces)

degree=zeros(face_number,1);
for i=1:face_number
    for j=1:3        
        degree(faces(i,j))=degree(faces(i,j))+1;                   
    end
end

