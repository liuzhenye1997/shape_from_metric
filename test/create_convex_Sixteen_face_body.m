%����һ��16���壬����������������м�8���������
%z�������������������㵽x-yƽ��ľ��룬����Ϊrand(1)*z+1
function points_sol=create_convex_Sixteen_face_body(z)
points_sol=zeros(10,3);

points_sol(1,1)=-(rand(1)+3);
points_sol(1,2)=0;
points_sol(2,1)=-(rand(1)+2);
points_sol(2,2)=(rand(1)+2);
points_sol(3,1)=0;
points_sol(3,2)=(rand(1)+3);
points_sol(4,1)=(rand(1)+2);
points_sol(4,2)=(rand(1)+2);
points_sol(5,1)=(rand(1)+3);
points_sol(5,2)=0;
points_sol(6,1)=(rand(1)+2);
points_sol(6,2)=-(rand(1)+2);
points_sol(7,1)=0;
points_sol(7,2)=-(rand(1)+3);
points_sol(8,1)=-(rand(1)+2);
points_sol(8,2)=-(rand(1)+2);

points_sol(9,3)=rand(1)*z+1;
points_sol(10,3)=-(rand(1)*z+1);

