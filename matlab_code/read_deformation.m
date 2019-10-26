function df=read_deformation(length,normal,face_number,face_hedge,points,vertex_theta,vertex_src,vertex_next,vertex_dst)
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
