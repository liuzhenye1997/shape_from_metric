%fix_-pi_pi_ambiguity
function vertex_A=fix_pi_pi_ambiguity(vertex_A,vertex_flip)
vertex_A_temp=vertex_A';
face_number=size(vertex_A,1);
for i=1:face_number
    for j=1:3
        if vertex_flip(3*i+j-3)~=1
            Aopp=vertex_A_temp(vertex_flip(3*(i-1)+j));
            if Aopp*vertex_A(i,j)>0 && 3*i+j-3<vertex_flip(3*i+j-3)
                vertex_A(i,j)=-vertex_A(i,j);
            end
        end
    end
end
