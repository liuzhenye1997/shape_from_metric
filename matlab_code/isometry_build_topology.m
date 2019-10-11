function [vertex_face,vertex_dst,vertex_next,vertex_prev,vertex_hedge,face_point,point_hedge]=isometry_build_topology(name)
[points,faces]=read_obj(name);
face_number=size(face,1);
point_number=size(points,1);
%create vertex attribute
vertex_face=[1:face_number;1:face_number;1:face_number]';
vertex_src=[faces(:,1),faces(:,3),faces(:,2)];
vertex_dst=[faces(:,3),faces(:,2),faces(:,1)];
vertex_next=[(0:face_number-1)*3+2;(0:face_number-1)*3+3;(0:face_number-1)*3+1]';
vertex_prev=[(0:face_number-1)*3+3;(0:face_number-1)*3+1;(0:face_number-1)*3+2]';
vertex_hedge=[(0:face_number-1)*3+1;(0:face_number-1)*3+2;(0:face_number-1)*3+3]';

%空缺半边结构及
%i@flip = hedge_srcvertex(0,
%             hedge_nextequiv(0,
%                 vertexhedge(0,@vtxnum)));
% if(@flip==@vtxnum){@flip=-1;} // boundary case

%create face attribute
face_hedge=3*(0:face_number-1)+1;
face_point=vertex_src(face_hedge);

%create point attribute
point_hedge=zeros(point_number);
for i=1:face_number
    for j=1:3
        point_hedge(face(i,j))=3*(i-1)+j;
    end
end