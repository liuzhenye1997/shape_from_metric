%This funvtion takes a triangle mesh with attributes (lambda) defined on faces describing triangle orientation quaternion, and a cached Laplace matrix , to produce a least-squares solution for the point positions.  This amounts to solving a Poisson problem. 
function points=isometry_possion_reconstruction(angle,vertex_prev,vertex_dst,vertex_src,vertex_flip,L_graph,length,lambda,vertex_face,vertex_theta)
point_number=size(L_graph,1);
face_number=size(lambda,1);
%vertex_edgevecΪ��İ�ߵ�����
vertex_edgevec=construction_edge_direction_per_face(lambda,vertex_face,vertex_theta);
vertex_edgevec_flip=vertex_edgevec(vertex_flip,:);
vertex_edgevec=0.5*(vertex_edgevec-vertex_edgevec_flip);
vertex_edgevec=vertex_edgevec.*length;
%points_div�ǵ��������㱾���ĺ͡������������������˽ǶȽ����˼�Ȩ��
points_div=compute_divergence(vertex_src,vertex_dst,point_number,face_number,vertex_edgevec,angle,vertex_prev);
%��õ��λ��
LII=L_graph(2:end,2:end);
points_divI=points_div(2:end,:);
solxyzI=LII\points_divI;
points(:,:)=[0 0 0;solxyzI];

%% �������������Ĳ�ļ�Ȩ��
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

%�������ѭ���滻���������һ���ִ���Ҳ�������С�
%��Ȼ����Ŀ���������ಢ���������ˡ�����������40��Ϊ�����������ʱ10s�࣬�������ʱ2s������
%����ԭ��̫���������������������Ҳ����Ҫʱ��ģ��ô��������̶�̫�͡�
%����iҲ����������ԭ����vertex_src�л�����ظ�Ԫ�أ���matlabֻ�����һ�β���������ʱ�޷���������
%��֮�ܶ�ط�������д�Ÿ���൫ȴû����д�ܿ�������Ϊ����ʱ�������
% points_div=zeros(point_number,3);
% for i=1:face_number
%     j=1:3;
%     points_div(vertex_src(i,:),:)=points_div(vertex_src(i,:),:)+weight(i,:)'.*vertex_edgevec(3*i+j-3,2:4);
%     points_div(vertex_dst(i,:),:)=points_div(vertex_dst(i,:),:)-weight(i,:)'.*vertex_edgevec(3*i+j-3,2:4);
% end
% points_divx=points_div(:,1);
% points_divy=points_div(:,2);
% points_divz=points_div(:,3);

%% ����ÿ����İ�ߵķ���
function vertex_edgevec=construction_edge_direction_per_face(lambda,vertex_face,vertex_theta)
qJ=[0 0 1 0];
vI=[1 0 0];
vertex_face_temp=vertex_face';
psi=lambda(vertex_face_temp(:),:);
vertex_theta_temp=vertex_theta';

result=Quaternion.batch_angleaxis_to_q(2.*vI.*vertex_theta_temp(:));
X=Quaternion.batch_multiplication(result,qJ);
vertex_edgevec=sandwich(psi,X,psi);

%% ����L�Ĺ������X����R��
function result=sandwich(L,X,R)
XR=Quaternion.batch_multiplication(X,R);
L(:,2:4)=-L(:,2:4);
result=Quaternion.batch_multiplication(L,XR);


