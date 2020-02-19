%输入x，即对角线长度和用来计算边长的原来的图形。
%输出error，即两条路径计算的点的位置的差和那个点的位置(result)。
function [error,result]=dis(x,points_sol)
faces=[1 2 5;2 3 5;1 4 5;3 4 5;1 2 6;2 3 6;1 4 6;3 4 6];
length=edge_length(faces,points_sol).^2;
l=zeros(size(points_sol,1),size(points_sol,1));
for i=1:size(faces,1)
    for j=1:3
        l(faces(i,mod(j-1,3)+1),faces(i,mod(j,3)+1))=length(i,j);
        l(faces(i,mod(j,3)+1),faces(i,mod(j-1,3)+1))=length(i,j);
    end
end
points=zeros(6,3);
points(5,:)=[0 0 x/2];
points(6,:)=[0 0 -x/2];

z1=(l(1,6)-l(1,5))*(2*x);
x1=-sqrt(l(1,5)-(z1-x/2)^2);
points(1,:)=[x1 0 z1];

z2=(l(2,6)-l(2,5))/(2*x);
k=(l(2,5)-(z2-x/2)^2)-(l(1,2)-(z2-z1)^2);
x2=(k+x1^2)/(2*x1);
y2=sqrt(l(2,5)-(z2-x/2)^2-x2^2);
points(2,:)=[x2 y2 z2];

z4=(l(4,6)-l(4,5))/(2*x);
k=(l(4,5)-(z4-x/2)^2)-(l(1,4)-(z4-z1)^2);
x4=(k+x1^2)/(2*x1);
y4=-sqrt(l(4,5)-(z4-x/2)^2-x4^2);
points(4,:)=[x4 y4 z4];

z3=(l(3,6)-l(3,5))/(2*x);
k=(l(3,5)-(z3-x/2)^2)-(l(2,3)-(z3-z2)^2);
c=l(3,5)-(z3-x/2)^2-((k+x2^2+y2^2)/(2*y2))^2;
b=x2*(k+x2^2+y2^2)/(y2^2);
a=1+(x2/y2)^2;
x3=(b+sqrt(b^2+4*a*c))/(2*a);
y3=(k+x2^2+y2^2-2*x2*x3)/(2*y2);
points(3,:)=[x3 y3 z3];

z3=(l(3,6)-l(3,5))/(2*x);
k=(l(3,5)-(z3-x/2)^2)-(l(4,3)-(z3-z4)^2);
c=l(3,5)-(z3-x/2)^2-((k+x4^2+y4^2)/(2*y4))^2;
b=x4*(k+x4^2+y4^2)/(y4^2);
a=1+(x4/y4)^2;
x3=(b+sqrt(b^2+4*a*c))/(2*a);
y3=(k+x2^2+y2^2-2*x2*x3)/(2*y2);
error=abs((x3-points(3,1)))^2+abs(y3-points(3,2))^2+abs(z3-points(3,3))^2;
result=[x3 y3 z3];
