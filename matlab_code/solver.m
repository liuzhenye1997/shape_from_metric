function [face_psi1re,face_psi1im,face_psi2re,face_psi2im,resx,resy,resz]=solver(mesh,iteration,vertex_theta,L_graph,angle,face_weight,face_src,face_dst,face_A,face_dihedral,face_theta,face_psi1re,face_psi1im,face_psi2re,face_psi2im)
%qexp_IA为τ=exp(i*A)表示的旋转转化成四元数的结果。
halfDihedralZ=0.5*face_dihedral(:).*[zeros(3*mesh.f_number,1) cos(face_theta(:)) sin(face_theta(:))];
face_A_temp=face_A';
IA=[face_A_temp(:) zeros(3*mesh.f_number,1) zeros(3*mesh.f_number,1)];
qexp_halfDihedralZ=Quaternion.batch_angleaxis_to_quaternion(-2*halfDihedralZ);
qexp_IA=Quaternion.batch_angleaxis_to_quaternion(-2*IA);
rho_inv=Quaternion.batch_quaternion_multiplication(qexp_halfDihedralZ,qexp_IA);
tau_inv=qexp_IA;
%z为论文中的z
face_theta_temp=face_theta';
z=vector2_expi(face_theta_temp(:));

face_src_temp=face_src';
face_dst_temp=face_dst';
psi=[face_psi1re(face_src_temp(:)) face_psi1im(face_src_temp(:)) face_psi2re(face_src_temp(:)) face_psi2im(face_src_temp(:))];
psi_dst=[face_psi1re(face_dst_temp(:)) face_psi1im(face_dst_temp(:)) face_psi2re(face_dst_temp(:)) face_psi2im(face_dst_temp(:))];
%mu为论文中的μ
mu=psi+Quaternion.batch_quaternion_multiplication(tau_inv,psi_dst);
mu_length=sqrt(sum(mu.^2,2));
avgfac=1./mu_length;
mu=mu./mu_length;
xi=Quaternion.batch_quaternion_multiplication(tau_inv,psi_dst)-psi;
%Z=z*j;I代表四元数i。
Z=[zeros(3*mesh.f_number,2) z];
I=[zeros(3*mesh.f_number,1) ones(3*mesh.f_number,1) zeros(3*mesh.f_number,2)];
I=Quaternion.batch_quaternion_multiplication(qexp_halfDihedralZ,I);
IZ=Quaternion.batch_quaternion_multiplication(I,Z);

Imu=Quaternion.batch_quaternion_multiplication(I,mu);
Zmu=Quaternion.batch_quaternion_multiplication(Z,mu);
IZmu=Quaternion.batch_quaternion_multiplication(IZ,mu);
F=tau_inv;
G=rho_inv;

%e_r,r_c,r_b,r_t代表论文中的μ，i*μ,z*j*μ,i*z*j*μ.
[k_real,k_connection,k_torsion,k_bending]=create_k(iteration,mesh.f_number);
e_r=mu;
e_c=Imu;
e_b=Zmu;
e_t=IZmu;
k_r=k_real;
k_c=k_connection;
k_t=k_torsion;
k_b=k_bending;

star1=sparse(1:3*mesh.f_number,1:3*mesh.f_number,face_weight);
M=sparse(1:mesh.f_number,1:mesh.f_number,mesh.f_area);
[L,M,DMU]=L_M_DMU(k_r,k_b,k_c,k_t,e_r,e_b,e_c,e_t,mesh.f_number,star1,mesh.f_ind,face_src,face_dst,F,G,M,avgfac);
%gradE_2nd_term即为算法5的g
gradE_2nd_term=grad_energy(DMU,xi,Z,I,IZ,mu,face_weight);
dt=Dt(iteration);
%利用梯度下降法最优化能量，即论文的Algorithm 5
[face_psi1re,face_psi1im,face_psi2re,face_psi2im]=gradient_flow(L,dt,mesh.f_number,M,gradE_2nd_term,face_psi1re,face_psi1im,face_psi2re,face_psi2im);
%由face_psi计算出点的位置
[resx,resy,resz,points]=isometry_possion_reconstruction(angle,mesh.v_prev,mesh.v_dst,mesh.v_src,mesh.v_flip,L_graph,mesh.length,face_psi2re,face_psi2im,face_psi1im,face_psi1re,mesh.v_face,vertex_theta);
%论文5.4.1,算法8.此处设置了a=0。
[face_psi1re,face_psi1im,face_psi2re,face_psi2im]=isometry_spinor_from_immersion(dt,mesh.length,mesh.f_number,vertex_theta,mesh.v_src,mesh.v_dst,mesh.v_next,mesh.f_hedge,points,mesh.faces,face_psi1re,face_psi1im,face_psi2re,face_psi2im);


%L，M为论文App D中的L和M。
function [L,M,DMU]=L_M_DMU(k_r,k_b,k_c,k_t,e_r,e_b,e_c,e_t,face_number,star1,face_ind,face_src,face_dst,F,G,M,avgfac)
w11 = k_r.*e_r(:,1).*e_r(:,1) + k_b.*e_b(:,1).*e_b(:,1) + k_c.*e_c(:,1).*e_c(:,1) + k_t.*e_t(:,1).*e_t(:,1);
w12 = k_r.*e_r(:,1).*e_r(:,2) + k_b.*e_b(:,1).*e_b(:,2) + k_c.*e_c(:,1).*e_c(:,2) + k_t.*e_t(:,1).*e_t(:,2);
w13 = k_r.*e_r(:,1).*e_r(:,3) + k_b.*e_b(:,1).*e_b(:,3) + k_c.*e_c(:,1).*e_c(:,3) + k_t.*e_t(:,1).*e_t(:,3);
w14 = k_r.*e_r(:,1).*e_r(:,4) + k_b.*e_b(:,1).*e_b(:,4) + k_c.*e_c(:,1).*e_c(:,4) + k_t.*e_t(:,1).*e_t(:,4);
w22 = k_r.*e_r(:,2).*e_r(:,2) + k_b.*e_b(:,2).*e_b(:,2) + k_c.*e_c(:,2).*e_c(:,2) + k_t.*e_t(:,2).*e_t(:,2);
w23 = k_r.*e_r(:,2).*e_r(:,3) + k_b.*e_b(:,2).*e_b(:,3) + k_c.*e_c(:,2).*e_c(:,3) + k_t.*e_t(:,2).*e_t(:,3);
w24 = k_r.*e_r(:,2).*e_r(:,4) + k_b.*e_b(:,2).*e_b(:,4) + k_c.*e_c(:,2).*e_c(:,4) + k_t.*e_t(:,2).*e_t(:,4);
w33 = k_r.*e_r(:,3).*e_r(:,3) + k_b.*e_b(:,3).*e_b(:,3) + k_c.*e_c(:,3).*e_c(:,3) + k_t.*e_t(:,3).*e_t(:,3);
w34 = k_r.*e_r(:,3).*e_r(:,4) + k_b.*e_b(:,3).*e_b(:,4) + k_c.*e_c(:,3).*e_c(:,4) + k_t.*e_t(:,3).*e_t(:,4);
w44 = k_r.*e_r(:,4).*e_r(:,4) + k_b.*e_b(:,4).*e_b(:,4) + k_c.*e_c(:,4).*e_c(:,4) + k_t.*e_t(:,4).*e_t(:,4);

W11=sparse(1:3*face_number,1:3*face_number,w11)*star1;
W12=sparse(1:3*face_number,1:3*face_number,w12)*star1;
W13=sparse(1:3*face_number,1:3*face_number,w13)*star1;
W14=sparse(1:3*face_number,1:3*face_number,w14)*star1;
W22=sparse(1:3*face_number,1:3*face_number,w22)*star1;
W23=sparse(1:3*face_number,1:3*face_number,w23)*star1;
W24=sparse(1:3*face_number,1:3*face_number,w24)*star1;
W33=sparse(1:3*face_number,1:3*face_number,w33)*star1;
W34=sparse(1:3*face_number,1:3*face_number,w34)*star1;
W44=sparse(1:3*face_number,1:3*face_number,w44)*star1;

W1=[W11,W12,W13,W14];
W2=[W12,W22,W23,W24];
W3=[W13,W23,W33,W34];
W4=[W14,W24,W34,W44];
W=[W1;W2;W3;W4];

face_src=face_src';
face_dst=face_dst';
d_row=[face_ind,face_ind]';
d_col=[face_src(:);face_dst(:)];

D1_val=[-ones(3*face_number,1),G(:,1)];
D2_val=[zeros(3*face_number,1),G(:,2)];
D3_val=[zeros(3*face_number,1),G(:,3)];
D4_val=[zeros(3*face_number,1),G(:,4)];


D1=sparse(d_row,d_col,D1_val,3*face_number,face_number);
D2=sparse(d_row,d_col,D2_val,3*face_number,face_number);
D3=sparse(d_row,d_col,D3_val,3*face_number,face_number);
D4=sparse(d_row,d_col,D4_val,3*face_number,face_number);

D_1=[D1,-D2,-D3,-D4];
D_2=[D2,D1,-D4,D3];
D_3=[D3,D4,D1,-D2];
D_4=[D4,-D3,D2,D1];
D=[D_1;D_2;D_3;D_4];
L=D'*W*D;


ZM=sparse(size(M,1),size(M,2));
MN_1=[M,ZM,ZM,ZM];
MN_2=[ZM,M,ZM,ZM];
MN_3=[ZM,ZM,M,ZM];
MN_4=[ZM,ZM,ZM,M];
M=[MN_1;MN_2;MN_3;MN_4];

AVG1_val=[avgfac;(F(:,1).*avgfac)];
AVG2_val=[zeros(3*face_number,1);(F(:,2).*avgfac)];
AVG3_val=[zeros(3*face_number,1);(F(:,3).*avgfac)];
AVG4_val=[zeros(3*face_number,1);(F(:,4).*avgfac)];

AVG1=sparse(d_row,d_col,AVG1_val,3*face_number,face_number);
AVG2=sparse(d_row,d_col,AVG2_val,3*face_number,face_number);
AVG3=sparse(d_row,d_col,AVG3_val,3*face_number,face_number);
AVG4=sparse(d_row,d_col,AVG4_val,3*face_number,face_number);
AVG_1=[AVG1,-AVG2,-AVG3,-AVG4];
AVG_2=[AVG2,AVG1,-AVG4,AVG3];
AVG_3=[AVG3,AVG4,AVG1,-AVG2];
AVG_4=[AVG4,-AVG3,AVG2,AVG1];
AVG=[AVG_1;AVG_2;AVG_3;AVG_4];

P11=sparse(1:3*face_number,1:3*face_number,1-e_r(:,1).*e_r(:,1));
P12=sparse(1:3*face_number,1:3*face_number,-e_r(:,1).*e_r(:,2));
P13=sparse(1:3*face_number,1:3*face_number,-e_r(:,1).*e_r(:,3));
P14=sparse(1:3*face_number,1:3*face_number,-e_r(:,1).*e_r(:,4));
P22=sparse(1:3*face_number,1:3*face_number,1-e_r(:,2).*e_r(:,2));
P23=sparse(1:3*face_number,1:3*face_number,-e_r(:,2).*e_r(:,3));
P24=sparse(1:3*face_number,1:3*face_number,-e_r(:,2).*e_r(:,4));
P33=sparse(1:3*face_number,1:3*face_number,1-e_r(:,3).*e_r(:,3));
P34=sparse(1:3*face_number,1:3*face_number,-e_r(:,3).*e_r(:,4));
P44=sparse(1:3*face_number,1:3*face_number,1-e_r(:,4).*e_r(:,4));
P_1=[P11,P12,P13,P14];
P_2=[P12,P22,P23,P24];
P_3=[P13,P23,P33,P34];
P_4=[P14,P24,P34,P44];
P=[P_1;P_2;P_3;P_4];

DMU=P*AVG;
%根据迭代次数调整k，即论文中的μ
function [k_real,k_connection,k_torsion,k_bending]=create_k(iteration,face_number)
if iteration<10
    k_real=ones(3*face_number,1);
    k_connection=ones(3*face_number,1);
    k_torsion=ones(3*face_number,1);
    k_bending=ones(3*face_number,1);
elseif iteration<20
    k_real=ones(3*face_number,1);
    k_connection=ones(3*face_number,1);
    k_torsion=ones(3*face_number,1)*0.1;
    k_bending=ones(3*face_number,1)*0.1;
elseif iteration<50
    k_real=ones(3*face_number,1);
    k_connection=ones(3*face_number,1);
    k_torsion=ones(3*face_number,1)*0.1;
    k_bending=ones(3*face_number,1)*0.01;
else
    k_real=ones(3*face_number,1);
    k_connection=ones(3*face_number,1);
    k_torsion=ones(3*face_number,1)*0.1;
    k_bending=zeros(3*face_number,1);
end
%根据迭代次数调整dt
function result=Dt(iteration)
if iteration<10
    result=2;
elseif iteration<20
    result=1;
elseif iteration<50
    result=1;
else
    result=0.2;
end
%利用梯度下降法最优化能量，即论文的Algorithm 5
function [points_psi1re,points_psi1im,points_psi2re,points_psi2im]=gradient_flow(L,dt,face_number,M,gradE_2nd_term,points_psi1re,points_psi1im,points_psi2re,points_psi2im)
%由于M是对角矩阵
area=M((size(M,1)+1)*(1:size(M,1))-size(M,1))';

Psi=[points_psi1re;points_psi1im;points_psi2re;points_psi2im];
Psi=Psi-dt.*(gradE_2nd_term./area);

LeftOp=M+dt.*L;
Psi=LeftOp\(M*Psi);

%归一化成单位向量
points_psi1re=Psi(1:face_number);
points_psi1im=Psi(face_number+1:2*face_number);
points_psi2re=Psi(2*face_number+1:3*face_number);
points_psi2im=Psi(3*face_number+1:4*face_number);
points_psi=[points_psi1re points_psi1im points_psi2re points_psi2im];
points_psi_norm=sqrt(sum(points_psi.^2,2));
points_psi1re=points_psi1re./points_psi_norm;
points_psi1im=points_psi1im./points_psi_norm;
points_psi2re=points_psi2re./points_psi_norm;
points_psi2im=points_psi2im./points_psi_norm;
%求gradE_2nd_term，即算法5的g
function gradE_2nd_term=grad_energy(DMU,xi,Z,I,IZ,mu,face_weight)
hp_1_xi=xi;
hp_Z_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(Z),xi);
hp_I_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(I),xi);
hp_IZ_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(IZ),xi);

ip_mu_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(mu),hp_1_xi);
ip_Zmu_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(mu),hp_Z_xi);
ip_Imu_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(mu),hp_I_xi);
ip_IZmu_xi=Quaternion.batch_quaternion_multiplication(Quaternion.qinvert(mu),hp_IZ_xi);

gradE_mu=face_weight.*(hp_1_xi.*ip_mu_xi(:,1)+hp_Z_xi.*ip_Zmu_xi(:,1)+hp_I_xi.*ip_Imu_xi(:,1)+hp_IZ_xi.*ip_IZmu_xi(:,1));
%gradE_2nd_term即为算法5中的g
gradE_mu_temp=[gradE_mu(:,1);gradE_mu(:,2);gradE_mu(:,3);gradE_mu(:,4)];
gradE_2nd_term=(DMU')*(gradE_mu_temp(:));

