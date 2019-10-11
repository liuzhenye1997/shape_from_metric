function [face_dihedral,face_weight,face_dual_length,face_length,face_theta,point_isOld,point_weight,face_dst,face_flip,face_src,face_A]=isometry_build_primal_graph(vertex_dual_length,vertex_src,vertex_flip,vertex_dst,point_area,point_number,length)
point_isOld=ones(size(point_number));
point_weight=point_area;

face_dst=vertex_dst;
face_flip=vertex_flip;
face_src=vertex_src;
face_A=0;
face_dihedral=0;
face_theta=0;
face_length=length;
face_dual_length=vertex_dual_length;
face_weight=vertex_dual_length./length;