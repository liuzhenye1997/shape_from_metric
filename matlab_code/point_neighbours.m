function nbrs=point_neighbours(face_number,faces,point_index)
nbrs=[];
for i=1:face_number
    for j=1:3
        if faces(i,j)==point_index
            nbrs=[nbrs faces(i,mod(j,3)+1) faces(i,mod(j+1,3)+1)];           
        end
    end
end

nbrs=sort(nbrs);
nbrs=nbrs(:,1:2:size(nbrs,2));
