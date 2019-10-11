function vertex_edgevec=construction_edge_direction_per_face(face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta)
qI=[0 1 0 0];
qJ=[0 0 1 0];
qK=[0 0 0 1];
vI=[1 0 0];
vertex_face_temp=vertex_face';
psi=[face_psi1re(vertex_face_temp(:)),face_psi1im(vertex_face_temp(:)),face_psi2re(vertex_face_temp(:)),face_psi2im(vertex_face_temp(:))];
vertex_theta_temp=vertex_theta';
result=batch_angleaxis_to_quaternion(2.*vI.*vertex_theta_temp(:));
X=batch_quaternion_multiplication(result,qJ);
qev=sandwich(psi,X,psi);
vertex_edgevec=qev;
