function vertex_isCotree=set_isCotree_for_opp_hedges(vertex_isCotree,face_number,vertex_flip)
for i=1:face_number
    for j=1:3
        if vertex_isCotree(3*i+j-3)==1
            vtx_flip=vertex_flip(3*i+j-3);
            if vtx_flip~=-1
                vertex_isCotree(vtx_flip)=-1;
            end
        end
    end
end