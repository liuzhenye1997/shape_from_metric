%计算laplacian矩阵。d0每一行代表一条边，-1代表出射点，1代表入射点。
%故d0'*d0即为一般的laplacian矩阵，这里使用了point_weight进行了加权。
function [L,M,d0,star1]=regular_laplacian(face_weight,point_weight,face_number,point_number,vertex_src,vertex_dst)
vertex_src=vertex_src';
vertex_dst=vertex_dst';
vertex_src=vertex_src(:);
vertex_dst=vertex_dst(:);
d0_row=[1:3*face_number,1:3*face_number];
d0_col=[vertex_src,vertex_dst];
d0_val=[-ones(face_number*3,1),ones(face_number*3,1)];
d0=sparse(d0_row,d0_col,d0_val,face_number*3,point_number);

star0_right=sparse(1:length(point_weight(:)),1:length(point_weight(:)),point_weight);
star1_right=sparse(1:3*face_number,1:3*face_number,face_weight);
L=-d0'*(star1_right*d0);
M=star0_right;
star1=star1_right;