%This node takes for every hedge a square root of the parallel transport of the tangent bundle.
%Since the parallel transport r = exp(i*vertex_A) is stored as an angle `vertex_A` per hedge, 
%this node effectively replaces `vertex_A` by `vertex_A/2` up to integer multiple of π, yielding a spin connection τ = exp(i*vertex_A) up to signs.  
%The signs are chosen so that τ defines a spin structure (all point is prescribed to be immersed).
%This is done with Algorithm 2 of Appendix A in the paper .
%通过调整τ的正负号使每个vertex star都是figure-0，其中τ对应论文中的τ。即论文算法2的第三四行，附录A
function [vertex_A,points_isBad]=isometry_bundle_sqrt(points_gauss_curvature,vertex_A,vertex_src,vertex_flip,point_number,faces)
face_number=size(faces,1);
points_holonomy=points_gauss_curvature;
points_holonomy=points_holonomy/2;
vertex_A=vertex_A/2;
points_holonomyFromA=zeros(point_number,1);

for i=1:face_number
    for j=1:3
        points_holonomyFromA(vertex_src(i,j))=points_holonomyFromA(vertex_src(i,j))+vertex_A(i,j);
    end
end

points_isBad=label_bad_vertices(points_holonomyFromA,points_holonomy,2);

points_local_isBdy=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        if vertex_flip(3*i+j-3)==-1
            points_local_isBdy(vertex_src(i,j))=1;
            points_local_isBdy(vertex_dst(i,j))=1;
        end
    end
end
root_candidate=find_root_candidate(points_local_isBdy);

coroot=1;
[~,~,~,vertex_isCotree]=face_breadth_first_search(vertex_flip,face_number,coroot);

vertex_isCotree=set_isCotree_for_opp_hedges(vertex_isCotree,face_number,vertex_flip);

root=root_candidate;
[points_distance,points_parent,~,~]=tree_breadth_first_search(point_number,faces,root,vertex_isCotree);

max_points_distance=max(points_distance);
points_cumBad=points_isBad;
for i=max_points_distance:-1:1
    d=i;
    for j=1:point_number
        if points_distance(j)==d
             points_cumBad(points_parent(j))=points_cumBad(points_parent(j))+points_cumBad(j);
        end
    end
end

degree=points_degree(face_number,faces);
max_degree=max(degree);
pointvertices=point_to_vertices(max_degree,point_number,faces);
vertex_flipCount=zeros(face_number*3,1);
vertex_src_temp=vertex_src';
for i=1:point_number
    vtx=1;
    vtx_opp=1;
    for j=1:degree(i)
        vtxtmp=pointvertices(i,j);
        vtxtmp_opp=vertex_flip(vtxtmp);
        if vertex_src_temp(vtxtmp_opp)==points_parent(i)
            vtx=vtxtmp;
            vtx_opp=vtxtmp_opp;
            break;
        end
    end
    vertex_flipCount(vtx)=points_cumBad(i);
    vertex_flipCount(vtx_opp)=points_cumBad(i);
end

for i=1:face_number
    for j=1:3
        if mod(vertex_flipCount(3*i+j-3),2)==1
            vertex_A(i,j)=vertex_A(i,j)+pi;
        end        
        vertex_A(i,j)=atan2(sin(vertex_A(i,j)),cos(vertex_A(i,j)));
    end
end

vertex_A_temp=vertex_A';
for i=1:face_number
    for j=1:3
        if vertex_flip(3*i+j-3)~=-1
            Aopp=vertex_A_temp(vertex_flip(3*i+j-3));
            if Aopp*vertex_A(i,j)>0 && 3*i+j-3<vertex_flip(3*i+j-3)
                vertex_A(i,j)=-1*vertex_A(i,j);
            end
        end
    end
end

points_holonomyFromA=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        points_holonomyFromA(vertex_src(i,j))=points_holonomyFromA(vertex_src(i,j))+vertex_A(i,j);
    end
end
evenOrodd=floor((points_holonomy(:)-points_holonomyFromA(:))/pi+0.5);
points_isBad=mod(evenOrodd,2);

%%
function vertex_isCotree=set_isCotree_for_opp_hedges(vertex_isCotree,face_number,vertex_flip)
for i=1:face_number
    for j=1:3
        if vertex_isCotree(3*i+j-3)==1
            vtx_flip=vertex_flip(3*i+j-3);
            if vtx_flip~=-1
                vertex_isCotree(vtx_flip)=-1;
            end
        end
    end
end

%%
function [points_distance,points_parent,points_parenting_vtx,vertex_isTree]=tree_breadth_first_search(point_number,faces,root,vertex_isCotree)
face_number=size(faces,1);
vertex_isTree=zeros(face_number*3,1);
points_distance=-ones(point_number,1);
points_parent=-ones(point_number,1);
points_parenting_vtx=-ones(point_number,1);
Q=[];
distance=-ones(point_number,1);
distance(root)=0;
points_distance(root)=0;
Q=[Q root];
degree=points_degree(face_number,faces);
max_degree=max(degree);

pointvertices=point_to_vertices(max_degree,point_number,faces);
pointvertices=pointvertices-1;

face_temp=faces';

while size(Q,2)>0
    curr=Q(1);
    Q=Q(2:size(Q,2));
    vtxs=pointvertices(curr,1:degree(curr));

    for i=1:degree(curr)
        j=mod(vtxs(i)+1,3)+1;
        k=floor(vtxs(i)/3+0.1);       
        nbr=face_temp(3*k+j);        
        if distance(nbr)==-1 && vertex_isCotree(vtxs(i)+1)==0
            distance(nbr)=distance(curr)+1;
            points_distance(nbr)=distance(nbr);
            points_parent(nbr)=curr;
            points_parenting_vtx(nbr)=vtxs(i)+1;
            vertex_isTree(vtxs(i)+1)=1;
            Q=[Q nbr];
        end        
    end
end

%%
function points_isBad=label_bad_vertices(points_holonomyFromA,points_holonomy,type)
point_number=size(points_holonomy,2);
points_isBad=zeros(point_number,1);
if type==1
    for i=1:point_number
        evenOrOdd=floor((points_holonomy(i)-points_holonomyFromA(i))/pi+0.5);
        points_isBad(i)=mod(evenOrOdd,2);
    end
else
    for i=1:point_number
        plusorminus=cos(points_holonomyFromA(i))*cos(points_holonomy(i))+sin(points_holonomyFromA(i))*sin(points_holonomy(i));
        points_isBad(i)=plusorminus<0;
    end
end

%% 由于vertex和point是多对一的关系。故pointvertices的每一行指的是point i所对应的vertex。
function pointvertices=point_to_vertices(max_degree,point_number,faces)
face_number=size(faces,1);
pointvertices=zeros(point_number,max_degree);
degree_temp=zeros(point_number,1);
for i=1:face_number
    for j=1:3
        degree_temp(faces(i,j))=degree_temp(faces(i,j))+1;
        pointvertices(faces(i,j),degree_temp(faces(i,j)))=3*i+j-3;    
    end
end

%%
function root_candidate=find_root_candidate(points_local_isBdy)
point_number=size(points_local_isBdy,1);
root_candidate=1;
for i=1:point_number
    isBdy=points_local_isBdy(i);
    if isBdy==1
        root_candidate=i;
        break;
    end
end

%%
function [face_distance,face_parent,face_parenting_vtx,vertex_isCotree]=face_breadth_first_search(vertex_flip,face_number,coroot)
Q=[];
face_distance=-ones(face_number,1);
face_parent=-ones(face_number,1);
face_parenting_vtx=-ones(face_number,1);
vertex_isCotree=zeros(3*face_number,1);
distance=-ones(face_number,1);
distance(coroot)=0;
Q=[Q coroot];
while size(Q,2)>0
    curr=Q(1);
    Q=Q(2:size(Q,2));
    for i=1:3
        vtx_flip=vertex_flip(3*curr+i-3);
        if vtx_flip==-1
            continue
        end
        nbr=floor(vtx_flip/3-0.1)+1;               
        if distance(nbr)==-1
            distance(nbr)=distance(curr)+1;
            face_distance(nbr)=distance(nbr);
            face_parent(nbr)=curr;
            face_parenting_vtx(nbr)=3*curr+i-3;
            vertex_isCotree(3*curr+i-3)=1;
            Q=[Q nbr];
        end
    end
end

