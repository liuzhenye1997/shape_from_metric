function [points_distance,points_parent]=breadth_first_search(faces,points_parent,points_distance,root)
Q=[];
face_number=size(faces,1);
point_number=size(points_distance,1);
distance=-ones(point_number,1);
distance(root)=0;
points_distance(root)=0;
Q=[Q root];
points_nbrs=points_neighbours(face_number,faces,point_number);


degree=points_degree(face_number,faces);

while size(Q,2)>0   
    curr=Q(1);
    Q=Q(2:size(Q,2));
    nbrs=points_nbrs(curr,:);
    for i=1:degree(curr)
        if distance(nbrs(i))==-1
            distance(nbrs(i))=distance(curr)+1;
            points_distance(nbrs(i))=distance(nbrs(i));
            points_parent(nbrs(i))=curr;
            Q=[Q nbrs(i)];
        end
    end
end