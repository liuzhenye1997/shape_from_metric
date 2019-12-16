function [solver1,lambda,L,index,M,P]=solver(solver1,g_is_zero,mesh,iteration,vertex_theta,L_graph,angle,face_weight,face_src,face_dst,face_A,face_dihedral,face_theta,lambda,L,index,M,P)
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

%Z=z*j;I代表四元数i。
Z=[zeros(3*mesh.f_number,2) z];
IZ=Quaternion.batch_multiplication(I,Z);

Imu=Quaternion.batch_multiplication(I,mu);
Zmu=Quaternion.batch_multiplication(Z,mu);
IZmu=Quaternion.batch_multiplication(IZ,mu);

%mu,Imu,Zmu,IZmu代表论文中的μ，i*μ,z*j*μ,i*z*j*μ.k则为论文中的ε
[k_r,k_c,k_t,k_b]=create_k(iteration);

% L，M为论文App D中的L和M。P和AVG用于求论文中的矩阵g
if iteration<1
    M=sparse(1:mesh.f_number*4,1:mesh.f_number*4,repelem(mesh.f_area,1,4));
end
[L,P,AVG,index]=L_P_AVG(iteration,g_is_zero,k_r,k_b,k_c,k_t,mu,Zmu,Imu,IZmu,mesh.f_number,face_weight,mesh.f_ind,face_dst,tau,rho,mesh.v_flip,avgfac,L,index,P);

%g即为算法5的g
if g_is_zero
    g=0;
else
    g=grad_energy(P,AVG,xi,Z,I,IZ,mu,Imu,Zmu,IZmu,face_weight);
end
dt=Dt(iteration);
%利用梯度下降法最优化能量，即论文的Algorithm 5
[solver1,lambda]=gradient_flow(iteration,solver1,L,dt,mesh.f_number,M,g,lambda);
%由lambda计算出点的位置
points=isometry_possion_reconstruction(angle,mesh.v_prev,mesh.v_dst,mesh.v_src,mesh.v_flip,L_graph,mesh.length,lambda,mesh.v_face,vertex_theta);
%论文5.4.1,算法8.此处设置了a=0。
lambda=isometry_spinor_from_immersion(dt,mesh.length,mesh.f_number,vertex_theta,mesh.v_src,mesh.v_dst,mesh.v_next,mesh.f_hedge,points,mesh.faces,lambda);
%画图
trimesh(mesh.faces, points(:,1), points(:,2), points(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
drawnow; pause(0.001);


%% L，M为论文App D中的L和M。P和AVG用于求论文中的矩阵g
function [L,P,AVG,index]=L_P_AVG(iteration,g_is_zero,k_r,k_b,k_c,k_t,mu,Zmu,Imu,IZmu,face_number,face_weight,face_ind,face_dst,tau,rho,face_flip,avgfac,L,index,P)
%w即为论文App D中的W
w=zeros(3*face_number,4,4);
for i=1:4
    for j=1:4
        w(:,i,j)=(k_r*mu(:,i).*mu(:,j) + k_b*Zmu(:,i).*Zmu(:,j) + k_c*Imu(:,i).*Imu(:,j) + k_t*IZmu(:,i).*IZmu(:,j)).*face_weight;
    end
end

face_dst_temp=face_dst';
face_dst_temp=face_dst_temp(:);
%D即为论文App D中的D
D=zeros(3*face_number,4,4);
D(:,1,:)=[rho(:,1) -rho(:,2) -rho(:,3) -rho(:,4)];
D(:,2,:)=[rho(:,2) rho(:,1) -rho(:,4) rho(:,3)];
D(:,3,:)=[rho(:,3) rho(:,4) rho(:,1) -rho(:,2)];
D(:,4,:)=[rho(:,4) -rho(:,3) rho(:,2) rho(:,1)];
%Q存储者L中的非零项的值。matrix_multiplication_c令两批4*4矩阵相乘
Q=zeros(3*face_number,4,4,4);
Q_temp=reshape(matrix_multiplication_c(D,w,true),3*face_number,4,4);
Q(:,1,:,:)=reshape(matrix_multiplication_c(Q_temp,D,false),3*face_number,4,4);
Q(:,2,:,:)=-Q_temp;
Q(:,3,:,:)=-reshape(matrix_multiplication_c(w,D,false),3*face_number,4,4);
Q(:,4,:,:)=w;
%L_row，L_col为L中非零项的行与列的下标
if iteration<1
    L_row=zeros(4,4,face_number,4);
    L_col=zeros(4,4,face_number,4);
    i=1:face_number;
    for j=1:4
        for k=1:4
            L_row(1,j,:,k)=i+(j-1)*face_number;
            L_col(1,k,:,j)=i+(j-1)*face_number;
            L_row(2,j,:,k)=face_dst(i,1)+(j-1)*face_number;
            L_col(2,k,:,j)=i+(j-1)*face_number;
            L_row(3,j,:,k)=face_dst(i,2)+(j-1)*face_number;
            L_col(3,k,:,j)=i+(j-1)*face_number;
            L_row(4,j,:,k)=face_dst(i,3)+(j-1)*face_number;
            L_col(4,k,:,j)=i+(j-1)*face_number;   
        end
    end
end
%将D中位置相同的非零项加起来，得到L_val
L_val=zeros(4,face_number,4,4);
face_flip=reshape(face_flip',3,face_number);
i=1:face_number ;
I=[3*i-2;3*i-1;3*i];
L_val(2:4,:,:,:)=reshape(Q(I,2,:,:)+Q(face_flip(:),3,:,:),3,face_number,4,4);
L_val(1,:,:,:)=Q(face_flip(1,:),1,:,:)+Q(face_flip(2,:),1,:,:)+Q(face_flip(3,:),1,:,:)+Q(3*i-2,4,:,:)+Q(3*i-1,4,:,:)+Q(3*i,4,:,:);
%为了与L_row，L_col的形式保持一致
L_val1=permute(L_val,[1,3,2,4]);
%论文App D中的矩阵L
if iteration<1
    L=sparse(L_row(:),L_col(:),L_val1(:)+10^-20,4*face_number,4*face_number);
    L_index=4*face_number.*(L_col(:)-1)+L_row(:);
    [~,index]=sort(L_index);
    index=int32(index)-1;
else
    update_matrix(L,L_val1+10^-20,index);
end


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
    if iteration<1
        spdiag=@(V)sparse(1:length(V(:)),1:length(V(:)),V(:));
        P11=spdiag(1-mu(:,1).*mu(:,1));
        P12=spdiag(-mu(:,1).*mu(:,2));
        P13=spdiag(-mu(:,1).*mu(:,3));
        P14=spdiag(-mu(:,1).*mu(:,4));
        P22=spdiag(1-mu(:,2).*mu(:,2));
        P23=spdiag(-mu(:,2).*mu(:,3));
        P24=spdiag(-mu(:,2).*mu(:,4));
        P33=spdiag(1-mu(:,3).*mu(:,3));
        P34=spdiag(-mu(:,3).*mu(:,4));
        P44=spdiag(1-mu(:,4).*mu(:,4));

        P=[P11,P12,P13,P14;...
           P12,P22,P23,P24;...
           P13,P23,P33,P34;...
           P14,P24,P34,P44];
    else
        P_val=[1-mu(:,1).*mu(:,1),-mu(:,1).*mu(:,2),-mu(:,1).*mu(:,3),-mu(:,1).*mu(:,4);...
                -mu(:,1).*mu(:,2),1-mu(:,2).*mu(:,2),-mu(:,2).*mu(:,3),-mu(:,2).*mu(:,4);...
                -mu(:,1).*mu(:,3),-mu(:,2).*mu(:,3),1-mu(:,3).*mu(:,3),-mu(:,3).*mu(:,4);...
                -mu(:,1).*mu(:,4),-mu(:,2).*mu(:,4),-mu(:,3).*mu(:,4),1-mu(:,4).*mu(:,4)]';    
        update_matrix(P,P_val);
    end
    %留作测试
%    spdiag=@(V)sparse(1:length(V(:)),1:length(V(:)),V(:));
%     P11=spdiag(1-mu(:,1).*mu(:,1));
%     P12=spdiag(-mu(:,1).*mu(:,2));
%     P13=spdiag(-mu(:,1).*mu(:,3));
%     P14=spdiag(-mu(:,1).*mu(:,4));
%     P22=spdiag(1-mu(:,2).*mu(:,2));
%     P23=spdiag(-mu(:,2).*mu(:,3));
%     P24=spdiag(-mu(:,2).*mu(:,4));
%     P33=spdiag(1-mu(:,3).*mu(:,3));
%     P34=spdiag(-mu(:,3).*mu(:,4));
%     P44=spdiag(1-mu(:,4).*mu(:,4));
%    
%     P=[P11,P12,P13,P14;...
%        P12,P22,P23,P24;...
%        P13,P23,P33,P34;...
%        P14,P24,P34,P44];
   %用时更长,约为上面的两倍。90%时间花在最后一行的构建稀疏矩阵上。
%     P_row=zeros(3*face_number,4,4);
%     P_col=zeros(3*face_number,4,4);
%     P_val=zeros(3*face_number,4,4);
%     for i=1:4
%         for j=1:4
%             P_row(:,i,j)=(1:3*face_number)+(i-1)*3*face_number;
%             P_col(:,i,j)=(1:3*face_number)+(j-1)*3*face_number;
%             if i==j
%                 P_val(:,i,j)=1-mu(:,i).*mu(:,j);
%             else
%                 P_val(:,i,j)=-mu(:,i).*mu(:,j);
%             end
%         end
%     end
%     P=sparse(P_row(:),P_col(:),P_val(:),3*face_number*4,3*face_number*4);
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
function [solver1,lambda]=gradient_flow(iteration,solver1,L,dt,face_number,M,g,lambda)
%由于M是对角矩阵
area=M((size(M,1)+1)*(1:size(M,1))-size(M,1))';

lambda=lambda(:);
lambda=lambda-dt.*(g./area);

LeftOp=M+dt.*L;
RightOp=M*lambda;
if iteration<1
    solver1 = splsolver(LeftOp, 'ldlt');
end
lambda=solver1.refactor_solve(LeftOp,RightOp);

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
gradE_mu_temp=gradE_mu(:);
g=AVG'*(P'*gradE_mu_temp);






