function [points_divx,points_divy,points_divz]=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev)
points_divx=zeros(point_number,1);
points_divy=zeros(point_number,1);
points_divz=zeros(point_number,1);
vertex_opp_angle=angle(vertex_prev);
weight=0.5./tan(vertex_opp_angle);
for i=1:face_number
    for j=1:3
        points_divx(vertex_src(i,j))=points_divx(vertex_src(i,j))+weight(i,j)*vertex_edgevec(3*i+j-3,2);
        points_divy(vertex_src(i,j))=points_divy(vertex_src(i,j))+weight(i,j)*vertex_edgevec(3*i+j-3,3);
        points_divz(vertex_src(i,j))=points_divz(vertex_src(i,j))+weight(i,j)*vertex_edgevec(3*i+j-3,4);
        points_divx(vertex_dst(i,j))=points_divx(vertex_dst(i,j))-weight(i,j)*vertex_edgevec(3*i+j-3,2);
        points_divy(vertex_dst(i,j))=points_divy(vertex_dst(i,j))-weight(i,j)*vertex_edgevec(3*i+j-3,3);
        points_divz(vertex_dst(i,j))=points_divz(vertex_dst(i,j))-weight(i,j)*vertex_edgevec(3*i+j-3,4);
    end
end
