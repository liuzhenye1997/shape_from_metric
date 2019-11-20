%This funvtion takes a triangle mesh with attributes (lambda) defined on faces describing triangle orientation quaternion, and a cached Laplace matrix , to produce a least-squares solution for the point positions.  This amounts to solving a Poisson problem. 
function points=isometry_possion_reconstruction(angle,vertex_prev,vertex_dst,vertex_src,vertex_flip,L_graph,length,lambda,vertex_face,vertex_theta)
point_number=size(L_graph,1);
face_number=size(lambda,1);
%vertex_edgevec为面的半边的向量
vertex_edgevec=construction_edge_direction_per_face(lambda,vertex_face,vertex_theta);
vertex_edgevec_flip=vertex_edgevec(vertex_flip,:);
vertex_edgevec=0.5*(vertex_edgevec-vertex_edgevec_flip);
vertex_edgevec=vertex_edgevec.*length;
%points_div是点的邻域与点本身差的和。在这样过程中利用了角度进行了加权。
points_div=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev);
%求得点的位置
LII=L_graph(2:end,2:end);
points_divI=points_div(2:end,:);
solxyzI=LII\points_divI;
points(:,:)=[0 0 0;solxyzI];

%% 计算点的邻域与点的差的加权和
function points_div=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev)
vertex_opp_angle=angle(vertex_prev);
weight=0.5./tan(vertex_opp_angle);
points_div=zeros(point_number,3);
for i=1:face_number
    for j=1:3
        points_div(vertex_src(i,j),1)=points_div(vertex_src(i,j),1)+weight(i,j)*vertex_edgevec(3*i+j-3,2);
        points_div(vertex_src(i,j),2)=points_div(vertex_src(i,j),2)+weight(i,j)*vertex_edgevec(3*i+j-3,3);
        points_div(vertex_src(i,j),3)=points_div(vertex_src(i,j),3)+weight(i,j)*vertex_edgevec(3*i+j-3,4);
        points_div(vertex_dst(i,j),1)=points_div(vertex_dst(i,j),1)-weight(i,j)*vertex_edgevec(3*i+j-3,2);
        points_div(vertex_dst(i,j),2)=points_div(vertex_dst(i,j),2)-weight(i,j)*vertex_edgevec(3*i+j-3,3);
        points_div(vertex_dst(i,j),3)=points_div(vertex_dst(i,j),3)-weight(i,j)*vertex_edgevec(3*i+j-3,4);
    end
end

%将上面的循环替换成下面的这一部分代码也可以运行。
%虽然下面的看起来更简洁并且向量化了。但是以运行40次为例，下面的用时10s多，上面的用时2s不到。
%具体原因不太清楚，可能是向量化本身也是需要时间的，该处向量化程度太低。
%不把i也给向量化的原因是vertex_src中会出现重复元素，而matlab只会进行一次操作。故暂时无法向量化。
%总之很多地方明显能写着更简洁但却没这样写很可能是因为所需时间更长。
% points_div=zeros(point_number,3);
% for i=1:face_number
%     j=1:3;
%     points_div(vertex_src(i,:),:)=points_div(vertex_src(i,:),:)+weight(i,:)'.*vertex_edgevec(3*i+j-3,2:4);
%     points_div(vertex_dst(i,:),:)=points_div(vertex_dst(i,:),:)-weight(i,:)'.*vertex_edgevec(3*i+j-3,2:4);
% end
% points_divx=points_div(:,1);
% points_divy=points_div(:,2);
% points_divz=points_div(:,3);

%% 计算每个面的半边的方向
function vertex_edgevec=construction_edge_direction_per_face(lambda,vertex_face,vertex_theta)
qJ=[0 0 1 0];
vI=[1 0 0];
vertex_face_temp=vertex_face';
psi=lambda(vertex_face_temp(:),:);
vertex_theta_temp=vertex_theta';

result=Quaternion.batch_angleaxis_to_q(2.*vI.*vertex_theta_temp(:));
X=Quaternion.batch_multiplication(result,qJ);
vertex_edgevec=sandwich(psi,X,psi);

%% 计算L的共轭乘以X乘以R。
function result=sandwich(L,X,R)
XR=Quaternion.batch_multiplication(X,R);
L(:,2:4)=-L(:,2:4);
result=Quaternion.batch_multiplication(L,XR);


