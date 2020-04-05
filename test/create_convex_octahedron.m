%创建一个16面体，由上下两个顶点和中间4个顶点组成
%z用来控制上下两个顶点到x-y平面的距离，距离为rand(1)*z+1
function points_sol=create_convex_octahedron(z)
points_sol=zeros(6,3);
points_sol(1,1)=-(rand(1)+1);
points_sol(1,3)=(rand(1)-0.5)*2;
points_sol(1,2)=(rand(1)-0.5)*2;
points_sol(2,2)=-(rand(1)+1);
points_sol(2,3)=(rand(1)-0.5)*2;
points_sol(2,1)=(rand(1)-0.5)*2;
points_sol(3,1)=rand(1)+1;
points_sol(3,3)=(rand(1)-0.5)*2;
points_sol(3,2)=(rand(1)-0.5)*2;
points_sol(4,2)=rand(1)+1;
points_sol(4,3)=(rand(1)-0.5)*2;
points_sol(4,1)=(rand(1)-0.5)*2;
points_sol(5,3)=rand(1)*z+1;
points_sol(6,3)=-(rand(1)*z+1);

