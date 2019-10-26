%This funvtion takes a triangle mesh with attributes (face_psi) defined on faces describing triangle orientation quaternion, and a cached Laplace matrix , to produce a least-squares solution for the point positions.  This amounts to solving a Poisson problem. 
function [resx,resy,resz,points]=isometry_possion_reconstruction(angle,vertex_prev,vertex_dst,vertex_src,vertex_flip,L_graph,length,face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta)
point_number=size(L_graph,1);
face_number=size(face_psi1im,1)*size(face_psi1im,2);
%vertex_edgevec为面的半边的向量
vertex_edgevec=construction_edge_direction_per_face(face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta);
vertex_edgevec_flip=vertex_edgevec(vertex_flip,:);
vertex_edgevec=0.5*(vertex_edgevec-vertex_edgevec_flip);
vertex_edgevec=vertex_edgevec.*length;
%points_div是点的邻域与点本身差的和。在这样过程中利用了角度进行了加权。
[points_divx,points_divy,points_divz]=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev);
%求得点的位置
LII=L_graph(2:end,2:end);
points_divxI=points_divx(2:end);
points_divyI=points_divy(2:end);
points_divzI=points_divz(2:end);
solxyzI=LII\[points_divxI points_divyI points_divzI];
solx=[0;solxyzI(:,1)];
soly=[0;solxyzI(:,2)];
solz=[0;solxyzI(:,3)];
points(:,:)=[solx,soly,solz];
%计算残差
resx=L_graph*solx-points_divx;
resy=L_graph*soly-points_divy;
resz=L_graph*solz-points_divz;


%计算点的邻域与点的差的加权和
function [points_divx,points_divy,points_divz]=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev)
points_divx=zeros(point_number,1);
points_divy=zeros(point_number,1);
points_divz=zeros(point_number,1);
vertex_opp_angle=angle(vertex_prev);
weight=0.5./tan(vertex_opp_angle);
for i=1:face_number
    for j=1:3
        points_divx(vertex_src(i,j))=points_divx(vertex_src(i,j))+weight(i,j)*vertex_edgevec(3*i+j-3,2);
        points_divy(vertex_src(i,j))=points_divy(vertex_src(i,j))+weight(i,j)*vertex_edgevec(3*i+j-3,3);
        points_divz(vertex_src(i,j))=points_divz(vertex_src(i,j))+weight(i,j)*vertex_edgevec(3*i+j-3,4);
        points_divx(vertex_dst(i,j))=points_divx(vertex_dst(i,j))-weight(i,j)*vertex_edgevec(3*i+j-3,2);
        points_divy(vertex_dst(i,j))=points_divy(vertex_dst(i,j))-weight(i,j)*vertex_edgevec(3*i+j-3,3);
        points_divz(vertex_dst(i,j))=points_divz(vertex_dst(i,j))-weight(i,j)*vertex_edgevec(3*i+j-3,4);
    end
end
%计算每个面的半边的方向
function vertex_edgevec=construction_edge_direction_per_face(face_psi2re,face_psi2im,face_psi1im,face_psi1re,vertex_face,vertex_theta)
qJ=[0 0 1 0];
vI=[1 0 0];
vertex_face_temp=vertex_face';
psi=[face_psi1re(vertex_face_temp(:)),face_psi1im(vertex_face_temp(:)),face_psi2re(vertex_face_temp(:)),face_psi2im(vertex_face_temp(:))];
vertex_theta_temp=vertex_theta';

result=Quaternion.batch_angleaxis_to_quaternion(2.*vI.*vertex_theta_temp(:));
X=Quaternion.batch_quaternion_multiplication(result,qJ);
vertex_edgevec=sandwich(psi,X,psi);
%计算L的共轭乘以X乘以R。
function result=sandwich(L,X,R)
XR=Quaternion.batch_quaternion_multiplication(X,R);
L(:,2:4)=-L(:,2:4);
result=Quaternion.batch_quaternion_multiplication(L,XR);


