%This node builds a graph with points corresponding to the faces of the input surface (in their `@ptnum` and `@primnum` indices) with primitives being polylines. 
%This node also carries over some hedge (Houdini-vertex) attributes from the original surface to the primitive attribute of the new geometry, including  `face_dihedral`, `face_theta`.
%This node also defines `weight` attribute on both points and primitives, as the point weight and edge weight for the graph, computed from the point area and the ratio of the primal length to the dual length stored on the original surface.  
function [face_dihedral,face_weight,face_vtx,face_dst,face_src]=isometry_build_dual_graph(lengths,face_number,vertex_flip,vertex_A,vertex_theta,vertex_dual_length)
%计算face_src,face_dst,face_vtx。其中face_src为顶点所在面，face_vtx为顶点的序号，face_dst为拥有同一条半边相邻的面。
face_src=zeros(face_number,3);
face_dst=zeros(face_number,3);
face_vtx=zeros(face_number,3);
for i=1:face_number
    for j=1:3
        flip=vertex_flip(3*i+j-3);
        %dsr_prim为第二维坐标
        dst_prim=(flip-(mod(flip-1,3)+1))/3+1;
        if flip>-1
            face_src(i,j)=i;
            face_dst(i,j)=dst_prim;
            face_vtx(i,j)=3*i+j-3;
        end
    end
end

% face_dihedral为面三条边与i-j平面的二面角。face_weight为面的三条半边的权重
face_dihedral=zeros(face_number,3);
face_weight=zeros(3*face_number,1);
for i=1:face_number
    for j=1:3
        length_temp=lengths(3*i+j-3);
        total_dual_length=vertex_dual_length(3*i+j-3)+vertex_dual_length(vertex_flip(3*i+j-3));
        face_weight(3*i+j-3)=0.5*length_temp/total_dual_length;
    end
end

