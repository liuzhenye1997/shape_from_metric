%计算vertex_flip,即顶点发出的半边的另一顶点。
function vertex_flip=vertex_flip(face_number,vertex_dst,vertex_src)
vertex_nexthedge=zeros(3*face_number,1);
vertex_dst_temp=vertex_dst';
vertex_src_temp=vertex_src';
vertex_dst_src=[vertex_dst_temp(:),vertex_src_temp(:)];
vertex_dst_src=sort(vertex_dst_src,2);
vertex_value=3*face_number*vertex_dst_src(:,1)+vertex_dst_src(:,2);
[s_vertex_value,index]=sort(vertex_value);
for i=1:3*face_number-1
    if s_vertex_value(i)==s_vertex_value(i+1)
        vertex_nexthedge(index(i))=index(i+1);
        vertex_nexthedge(index(i+1))=index(i);
    end
end

vertex_flip=vertex_nexthedge;