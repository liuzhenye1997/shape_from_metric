function [face_distance,face_parent,face_parenting_vtx,vertex_isCotree]=face_breadth_first_search(vertex_flip,face_number,coroot)
Q=[];
face_distance=-ones(face_number,1);
face_parent=-ones(face_number,1);
face_parenting_vtx=-ones(face_number,1);
vertex_isCotree=zeros(3*face_number,1);
distance=-ones(face_number,1);
distance(coroot)=0;
Q=[Q coroot];
while size(Q,2)>0
    curr=Q(1);
    Q=Q(2:size(Q,2));
    for i=1:3
        vtx_flip=vertex_flip(3*curr+i-3);
        if vtx_flip==-1
            continue
        end
        nbr=floor(vtx_flip/3-0.1)+1;
                
        if distance(nbr)==-1
            distance(nbr)=distance(curr)+1;
            face_distance(nbr)=distance(nbr);
            face_parent(nbr)=curr;
            face_parenting_vtx(nbr)=3*curr+i-3;
            vertex_isCotree(3*curr+i-3)=1;
            Q=[Q nbr];
        end
    end
end