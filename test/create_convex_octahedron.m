function points_sol=create_convex_octahedron()
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
points_sol(5,3)=rand(1)+1;
points_sol(6,3)=-(rand(1)+1);

