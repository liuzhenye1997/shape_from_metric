%论文5.4.1,算法8.此处设置了a=0。
function [face_psi1re,face_psi1im,face_psi2re,face_psi2im]=isometry_spinor_from_immersion(dt,length,face_number,vertex_theta,vertex_src,vertex_dst,vertex_next,face_hedge,points,faces,face_psi1re,face_psi1im,face_psi2re,face_psi2im)
normal=compute_normal(points,faces);

hedge0=face_hedge';
vertex_next_temp=vertex_next';
hedge1=vertex_next_temp(hedge0);
vertex_src_temp=vertex_src';
vertex_dst_temp=vertex_dst';
P_hedge0src=points(vertex_src_temp(hedge0),:);
P_hedge0dst=points(vertex_dst_temp(hedge0),:);
P_hedge1src=points(vertex_src_temp(hedge1),:);
P_hedge1dst=points(vertex_dst_temp(hedge1),:);
V0=P_hedge0dst-P_hedge0src;
V1=P_hedge1dst-P_hedge1src;

l0=length(hedge0);
l1=length(hedge1);
length_V0=sqrt(sum(V0.^2,2));
N=normal;
E=V0./length_V0;
NE=[N(:,2).*E(:,3)-N(:,3).*E(:,2) N(:,3).*E(:,1)-N(:,1).*E(:,3) N(:,1).*E(:,2)-N(:,2).*E(:,1)];
u=sum(E.*V1,2);
v=sum(NE.*V1,2);
vertex_theta_temp=vertex_theta';
theta=vertex_theta_temp(hedge1);
a=(u-l1.*l0./length_V0.*cos(theta))./(l1.*sin(theta));
b=v./(l1.*sin(theta));
F=a.*E+b.*NE;
df=zeros(face_number,9);
df(:,1:3)=N;
df(:,4:6)=length_V0./l0.*E;
df(:,7:9)=F;

q=zeros(face_number,4);
for i=1:face_number
    df_temp=[df(i,1:3);df(i,4:6);df(i,7:9)];
    P=sqrtm(df_temp'*df_temp);
    U=df_temp/P;
    q(i,:)=Quaternion.matrix_to_quaternion(U);
end
q(:,2:4)=-q(:,2:4);

face_psi1re_pre=face_psi1re;
face_psi1im_pre=face_psi1im;
face_psi2re_pre=face_psi2re;
face_psi2im_pre=face_psi2im;

prod=face_psi1re_pre.*q(:,1)+face_psi1im_pre.*q(:,2)+face_psi2re_pre.*q(:,3)+face_psi2im_pre.*q(:,4);
for i=1:face_number
    if prod(i)<0
        q(i,:)=-q(i,:);        
    end
end

face_psi1re_pre=face_psi1re;
face_psi2re_pre=face_psi2re;
face_psi2im_pre=face_psi2im;
face_psi1im_pre=face_psi1im;

face_psi1re=q(:,1);
face_psi1im=q(:,2);
face_psi2re=q(:,3);
face_psi2im=q(:,4);

%由于a=0，故事实上下面并没有改变face_psi;
a=0;
[face_psi2im,face_psi2re,face_psi1re,face_psi1im]=penalty_parameter(face_psi1re_pre,face_psi2re_pre,face_psi2im_pre,face_psi1im_pre,face_psi2im,face_psi2re,face_psi1re,face_psi1im,a,dt);


%算法8第7,8行
function [face_psi2im,face_psi2re,face_psi1re,face_psi1im]=penalty_parameter(face_psi1re_pre,face_psi2re_pre,face_psi2im_pre,face_psi1im_pre,face_psi2im,face_psi2re,face_psi1re,face_psi1im,a,dt)
if abs(a)<1*exp(-6)
    r=0;
else
    r=exp(-dt/a);
end

face_psi1im=face_psi1im*(1-r);
face_psi1re=face_psi1re*(1-r);
face_psi2re=face_psi2re*(1-r);
face_psi2im=face_psi2im*(1-r);

face_psi1im=face_psi1im+r*face_psi1im_pre;
face_psi2im=face_psi2im+r*face_psi2im_pre;
face_psi2re=face_psi2re+r*face_psi2re_pre;
face_psi1re=face_psi1re+r*face_psi1re_pre;

face_psi=[face_psi1im face_psi1re face_psi2im face_psi2re];
face_psi_norm=sqrt(sum(face_psi.^2,2));
face_psi1im=face_psi1im./face_psi_norm;
face_psi1re=face_psi1re./face_psi_norm;
face_psi2im=face_psi2im./face_psi_norm;
face_psi2re=face_psi2re./face_psi_norm;