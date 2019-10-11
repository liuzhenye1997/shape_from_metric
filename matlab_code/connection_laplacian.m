function [L,M,d0,star1]=connection_laplacian(face_weight,point_weight,face_number,point_number,vertex_src,vertex_dst,vertex_A)
vertex_A=vertex_A';
vertex_A=vertex_A(:);
d0_row=[1:3*face_number,1:3*face_number];
vertex_src=vertex_src';
vertex_dst=vertex_dst';
d0_col=[vertex_src(:);vertex_dst(:)];
d0_val=[-ones(face_number*3,1);exp(-1i*vertex_A)];
d0=sparse(d0_row,d0_col,d0_val,face_number*3,point_number);

star0=sparse(1:face_number,1:face_number,point_weight);
star1=sparse(1:3*face_number,1:3*face_number,face_weight);
L=-(d0')*star1*d0;
M=star0;