function lambda=solver(g_is_zero,mesh,iteration,vertex_theta,L_graph,angle,face_weight,face_src,face_dst,face_A,face_dihedral,face_theta,lambda)
%τ=exp(i*A)表示的旋转转化成四元数的结果。
halfDihedralZ=0.5*face_dihedral(:).*[zeros(3*mesh.f_number,1) cos(face_theta(:)) sin(face_theta(:))];
face_A_temp=face_A';
IA=[face_A_temp(:) zeros(3*mesh.f_number,2)];
I=[zeros(3*mesh.f_number,1) ones(3*mesh.f_number,1) zeros(3*mesh.f_number,2)];
qexp_halfDihedralZ=Quaternion.batch_angleaxis_to_q(-2*halfDihedralZ);
tau=Quaternion.batch_angleaxis_to_q(-2*IA);
rho=Quaternion.batch_multiplication(qexp_halfDihedralZ,tau);

%z为论文(14)式中的z
face_theta_temp=face_theta';
z=vector2_expi(face_theta_temp(:));

face_src_temp=face_src';
face_dst_temp=face_dst';
lambda_src=lambda(face_src_temp(:),:);
lambda_dst=lambda(face_dst_temp(:),:);
%mu为论文中的μ
mu=lambda_src+Quaternion.batch_multiplication(tau,lambda_dst);
mu_length=sqrt(sum(mu.^2,2));
avgfac=1./mu_length;
mu=mu./mu_length;
xi=Quaternion.batch_multiplication(tau,lambda_dst)-lambda_src;
% sum(abs(xi(:)))
%Z=z*j;I代表四元数i。
Z=[zeros(3*mesh.f_number,2) z];
IZ=Quaternion.batch_multiplication(I,Z);

Imu=Quaternion.batch_multiplication(I,mu);
Zmu=Quaternion.batch_multiplication(Z,mu);
IZmu=Quaternion.batch_multiplication(IZ,mu);

%mu,Imu,Zmu,IZmu代表论文中的μ，i*μ,z*j*μ,i*z*j*μ.k则为论文中的ε
[k_r,k_c,k_t,k_b]=create_k(iteration);
% L，M为论文App D中的L和M。P和AVG用于求论文中的矩阵g
[L,M,P,AVG]=L_M_P_AVG(g_is_zero,k_r,k_b,k_c,k_t,mu,Zmu,Imu,IZmu,mesh.f_number,face_weight,mesh.f_ind,face_src,face_dst,tau,rho,mesh.f_area,avgfac);

%g即为算法5的g
if g_is_zero
    g=0;
else
    g=grad_energy(P,AVG,xi,Z,I,IZ,mu,Imu,Zmu,IZmu,face_weight);
end
dt=Dt(iteration);
%利用梯度下降法最优化能量，即论文的Algorithm 5
lambda=gradient_flow(L,dt,mesh.f_number,M,g,lambda);
%由lambda计算出点的位置
points=isometry_possion_reconstruction(angle,mesh.v_prev,mesh.v_dst,mesh.v_src,mesh.v_flip,L_graph,mesh.length,lambda,mesh.v_face,vertex_theta);
%论文5.4.1,算法8.此处设置了a=0。
lambda=isometry_spinor_from_immersion(dt,mesh.length,mesh.f_number,vertex_theta,mesh.v_src,mesh.v_dst,mesh.v_next,mesh.f_hedge,points,mesh.faces,lambda);
%画图
trimesh(mesh.faces, points(:,1), points(:,2), points(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
drawnow; pause(0.001);


%% L，M为论文App D中的L和M。P和AVG用于求论文中的矩阵g
function [L,M,P,AVG]=L_M_P_AVG(g_is_zero,k_r,k_b,k_c,k_t,mu,Zmu,Imu,IZmu,face_number,face_weight,face_ind,face_src,face_dst,tau,rho,face_area,avgfac)
w=zeros(4,4,3*face_number);
for i=1:4
    for j=i:4
        w(i,j,:)=(k_r*mu(:,i).*mu(:,j) + k_b*Zmu(:,i).*Zmu(:,j) + k_c*Imu(:,i).*Imu(:,j) + k_t*IZmu(:,i).*IZmu(:,j)).*face_weight;
    end
end

spdiag=@(N,V)sparse(1:N,1:N,V(:));
W11=spdiag(3*face_number,w(1,1,:));
W12=spdiag(3*face_number,w(1,2,:));
W13=spdiag(3*face_number,w(1,3,:));
W14=spdiag(3*face_number,w(1,4,:));
W22=spdiag(3*face_number,w(2,2,:));
W23=spdiag(3*face_number,w(2,3,:));
W24=spdiag(3*face_number,w(2,4,:));
W33=spdiag(3*face_number,w(3,3,:));
W34=spdiag(3*face_number,w(3,4,:));
W44=spdiag(3*face_number,w(4,4,:));

face_src_temp=face_src';
face_dst_temp=face_dst';
d_row=[face_ind,face_ind]';
d_col=[face_src_temp(:);face_dst_temp(:)];

D1_val=[-ones(3*face_number,1),rho(:,1)];
D2_val=[zeros(3*face_number,1),rho(:,2)];

D1=sparse(d_row,d_col,D1_val,3*face_number,face_number);
D2=sparse(d_row,d_col,D2_val,3*face_number,face_number);

W1=[W11,W12;...
    W12,W22];
W2=[W13,W23;...
    W14,W24];
W3=[W33,W34;...
    W34,W44];
D11=[D1,-D2;...
    D2,D1];
L1=D11'*W1*D11;
L2=D11'*W2*D11;
L3=D11'*W3*D11;
L=[L1,L2';...
   L2,L3];

M=sparse(1:face_number*4,1:face_number*4,repelem(face_area,1,4));
if ~g_is_zero
    AVG1=sparse(face_ind,face_dst_temp(:),tau(:,1).*avgfac,3*face_number,face_number);
    AVG2=sparse(face_ind,face_dst_temp(:),tau(:,2).*avgfac,3*face_number,face_number);
    AVG3=sparse(face_ind,face_dst_temp(:),tau(:,3).*avgfac,3*face_number,face_number);
    AVG4=sparse(face_ind,face_dst_temp(:),tau(:,4).*avgfac,3*face_number,face_number);

    AVG=[AVG1,-AVG2,-AVG3,-AVG4;...
        AVG2,AVG1,-AVG4,AVG3;...
        AVG4,-AVG3,AVG2,AVG1;...
        AVG4,-AVG3,AVG2,AVG1];

    %P为论文App D中g的一组成部分
    P11=spdiag(3*face_number,1-mu(:,1).*mu(:,1));
    P12=spdiag(3*face_number,-mu(:,1).*mu(:,2));
    P13=spdiag(3*face_number,-mu(:,1).*mu(:,3));
    P14=spdiag(3*face_number,-mu(:,1).*mu(:,4));
    P22=spdiag(3*face_number,1-mu(:,2).*mu(:,2));
    P23=spdiag(3*face_number,-mu(:,2).*mu(:,3));
    P24=spdiag(3*face_number,-mu(:,2).*mu(:,4));
    P33=spdiag(3*face_number,1-mu(:,3).*mu(:,3));
    P34=spdiag(3*face_number,-mu(:,3).*mu(:,4));
    P44=spdiag(3*face_number,1-mu(:,4).*mu(:,4));

    P=[P11,P12,P13,P14;...
       P12,P22,P23,P24;...
       P13,P23,P33,P34;...
       P14,P24,P34,P44];
else 
    P=0;
    AVG=0;
end

%% 根据迭代次数调整k，即论文中的ε
function [k_real,k_connection,k_torsion,k_bending]=create_k(iteration)
if iteration<10
    k_real=1;
    k_connection=1;
    k_torsion=1;
    k_bending=1;
elseif iteration<20
    k_real=1;
    k_connection=1;
    k_torsion=0.1;
    k_bending=0.1;
elseif iteration<50
    k_real=1;
    k_connection=1;
    k_torsion=0.1;
    k_bending=0.01;
else
    k_real=1;
    k_connection=1;
    k_torsion=1;
    k_bending=0;
end

%% 根据迭代次数调整dt
function dt=Dt(iteration)
if iteration<10
    dt=2;
elseif iteration<20
    dt=1;
elseif iteration<50
    dt=1;
elseif iteration<70
    dt=0.2;
elseif iteration<100
    dt=0.1;   
elseif iteration
    dt=10^-1*10^-((iteration-100)/50);     
end

%% 利用梯度下降法最优化能量，即论文的Algorithm 5
function lambda=gradient_flow(L,dt,face_number,M,g,lambda)
%由于M是对角矩阵
area=M((size(M,1)+1)*(1:size(M,1))-size(M,1))';

lambda=lambda(:);
lambda=lambda-dt.*(g./area);

LeftOp=M+dt.*L;
RightOp=M*lambda;
solver = splsolver(LeftOp, 'ldlt');
lambda = solver.refactor_solve(LeftOp,RightOp);

%归一化成单位向量
lambda=[lambda(1:face_number) lambda(face_number+1:2*face_number) lambda(2*face_number+1:3*face_number) lambda(3*face_number+1:4*face_number)];
lambda_norm=sqrt(sum(lambda.^2,2));
lambda=lambda./lambda_norm;

%% 求gradE_2nd_term，即算法5的g
function g=grad_energy(P,AVG,xi,Z,I,IZ,mu,Imu,Zmu,IZmu,face_weight)
mu_xi=sum(mu.*xi,2);
Zmu_xi=sum(Zmu.*xi,2);
Imu_xi=sum(Imu.*xi,2);
IZmu_xi=sum(IZmu.*xi,2);
Z_xi=Quaternion.batch_multiplication(Quaternion.qinvert(Z),xi);
I_xi=Quaternion.batch_multiplication(Quaternion.qinvert(I),xi);
IZ_xi=Quaternion.batch_multiplication(Quaternion.qinvert(IZ),xi);

gradE_mu=face_weight.*(xi.*mu_xi+Z_xi.*Zmu_xi+I_xi.*Imu_xi+IZ_xi.*IZmu_xi);
%g即为算法5中的g
gradE_mu_temp=[gradE_mu(:,1);gradE_mu(:,2);gradE_mu(:,3);gradE_mu(:,4)];
g=AVG'*(P'*gradE_mu_temp);


