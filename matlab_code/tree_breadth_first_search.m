function [points_distance,points_parent,points_parenting_vtx,vertex_isTree]=tree_breadth_first_search(point_number,faces,root,vertex_isCotree)

faces=[faces(:,1),faces(:,3),faces(:,2)];
faces(1,:)

face_number=size(faces,1);
vertex_isTree=zeros(face_number*3,1);
points_distance=-ones(point_number,1);
points_parent=-ones(point_number,1);
points_parenting_vtx=-ones(point_number,1);
Q=[];
distance=-ones(point_number,1);
distance(root)=0;
points_distance(root)=0;
Q=[Q root];
degree=points_degree(face_number,faces);
max_degree=max(degree);

pointvertices=point_to_vertices(max_degree,point_number,faces);
pointvertices=pointvertices-1;
% pointvertices=zeros(point_number,max_degree);
% degree_temp=zeros(point_number,1);
% for i=1:face_number
%     for j=1:3
%         degree_temp(faces(i,j))=degree_temp(faces(i,j))+1;
%         pointvertices(faces(i,j),degree_temp(faces(i,j)))=3*i+j-4;    
%     end
% end

face_temp=faces';
face_temp(:,1)


while size(Q,2)>0
    curr=Q(1);
    Q=Q(2:size(Q,2));
    vtxs=pointvertices(curr,1:degree(curr));

    for i=1:degree(curr)
        j=mod(vtxs(i)+1,3)+1;
        k=floor(vtxs(i)/3+0.1);   
    
        nbr=face_temp(3*k+j);
        
        if distance(nbr)==-1 && vertex_isCotree(vtxs(i)+1)==0
            distance(nbr)=distance(curr)+1;
            points_distance(nbr)=distance(nbr);
            points_parent(nbr)=curr;
            points_parenting_vtx(nbr)=vtxs(i)+1;
            vertex_isTree(vtxs(i)+1)=1;
            Q=[Q nbr];
        end        
    end
end
