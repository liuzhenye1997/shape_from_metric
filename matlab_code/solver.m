function [face_psi1re,face_psi1im,face_psi2re,face_psi2im,resx,resy,resz,points_old]=solver(iteration,face_hedge,vertex_theta,vertex_face,L3,vertex_flip,vertex_src,vertex_dst,vertex_next,vertex_prev,angle,faces,points_psi1re,points_psi1im,points_psi2re,points_psi2im,face_weight,M,face_ind,star1,face_src,face_dst,face_A,length,points,face_number,face_dihedral,face_theta,face_psi1re,face_psi1im,face_psi2re,face_psi2im)
[Prev_k_real,Prev_k_connection,Prev_k_torsion,Prev_k_bending]=create_k(iteration,face_number);


halfDihedralZ=0.5*face_dihedral(:).*[zeros(3*face_number,1) cos(face_theta(:)) sin(face_theta(:))];
face_A_temp=face_A';
IA=[face_A_temp(:) zeros(3*face_number,1) zeros(3*face_number,1)];
qexp_halfDihedralZ=batch_angleaxis_to_quaternion(-2*halfDihedralZ);
qexp_IA=batch_angleaxis_to_quaternion(-2*IA);
rho_inv=batch_quaternion_multiplication(qexp_halfDihedralZ,qexp_IA);
tau_inv=qexp_IA;

face_theta_temp=face_theta';
z=vector2_expi(face_theta_temp(:));

face_src_temp=face_src';
face_dst_temp=face_dst';
psi=[face_psi1re(face_src_temp(:)) face_psi1im(face_src_temp(:)) face_psi2re(face_src_temp(:)) face_psi2im(face_src_temp(:))];
psi_dst=[face_psi1re(face_dst_temp(:)) face_psi1im(face_dst_temp(:)) face_psi2re(face_dst_temp(:)) face_psi2im(face_dst_temp(:))];

mu=psi+batch_quaternion_multiplication(tau_inv,psi_dst);

mu_length=sqrt(sum(mu.^2,2));
avgfac=1./mu_length;
mu=mu./mu_length;
xi=batch_quaternion_multiplication(tau_inv,psi_dst)-psi;

Z=[zeros(3*face_number,2) z];
I=[zeros(3*face_number,1) ones(3*face_number,1) zeros(3*face_number,2)];
I=batch_quaternion_multiplication(qexp_halfDihedralZ,I);
IZ=batch_quaternion_multiplication(I,Z);

Imu=batch_quaternion_multiplication(I,mu);
Zmu=batch_quaternion_multiplication(Z,mu);
IZmu=batch_quaternion_multiplication(IZ,mu);
F=tau_inv;
G=rho_inv;

e_r=mu;
e_c=Imu;
e_b=Zmu;
e_t=IZmu;
k_r=Prev_k_real;
k_c=Prev_k_connection;
k_t=Prev_k_torsion;
k_b=Prev_k_bending;


[L,MM,DMU]=L_MM_DMU(k_r,k_b,k_c,k_t,e_r,e_b,e_c,e_t,face_number,star1,face_ind,face_src,face_dst,F,G,M,avgfac);

gradE_mu=grad_energy_wrt_mu(xi,Z,I,IZ,mu,face_weight);
dt=Dt(iteration);
M=MM;
[points_psi1re,points_psi1im,points_psi2re,points_psi2im]=gradient_flow(L,dt,face_number,M,DMU,gradE_mu,points_psi1re,points_psi1im,points_psi2re,points_psi2im);
points_psi=[points_psi1re points_psi1im points_psi2re points_psi2im];
points_psi_norm=sqrt(sum(points_psi.^2,2));

points_psi1re=points_psi1re./points_psi_norm;
points_psi1im=points_psi1im./points_psi_norm;
points_psi2re=points_psi2re./points_psi_norm;
points_psi2im=points_psi2im./points_psi_norm;

face_psi1im=points_psi1im;
face_psi1re=points_psi1re;
face_psi2im=points_psi2im;
face_psi2re=points_psi2re;

[resx,resy,resz,points_old,points]=isometry_possion_reconstruction(faces,angle,vertex_prev,vertex_dst,vertex_src,vertex_flip,points,L3,length,face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta);

face_psi1re_pre=face_psi1re;
face_psi1im_pre=face_psi1im;
face_psi2re_pre=face_psi2re;
face_psi2im_pre=face_psi2im;
[face_psi1re,face_psi1im,face_psi2re,face_psi2im]=isometry_spinor_from_immersion(length,face_number,vertex_theta,vertex_src,vertex_dst,vertex_next,face_hedge,points,faces,face_psi1re,face_psi1im,face_psi2re,face_psi2im);

a=0;
[face_psi2im,face_psi2re,face_psi1re,face_psi1im]=penalty_parameter(face_psi1re_pre,face_psi2re_pre,face_psi2im_pre,face_psi1im_pre,face_psi2im,face_psi2re,face_psi1re,face_psi1im,a,dt);


