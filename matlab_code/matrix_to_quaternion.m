function q=matrix_to_quaternion(m)
w_temp=sqrt(m(1,1)+m(2,2)+m(3,3)+1)/2;
x_temp=sqrt(m(1,1)-m(2,2)-m(3,3)+1)/2;
y_temp=sqrt(-m(1,1)+m(2,2)-m(3,3)+1)/2;
z_temp=sqrt(-m(1,1)-m(2,2)+m(3,3)+1)/2;
MAX=max([w_temp x_temp y_temp z_temp]);
if w_temp==MAX
    w=w_temp;
    x=(m(2,3)-m(3,2))/(4*w);
    y=(m(3,1)-m(1,3))/(4*w);
    z=(m(1,2)-m(2,1))/(4*w);
elseif x_temp==MAX
    x=x_temp;
    w=(m(2,3)-m(3,2))/(4*x);
    y=(m(1,2)+m(2,1))/(4*x);
    z=(m(3,1)+m(1,3))/(4*x);
elseif y_temp==MAX
    y=y_temp;
    w=(m(3,1)-m(1,3))/(4*y);
    x=(m(1,2)+m(2,1))/(4*y);
    z=(m(2,3)+m(3,2))/(4*y);
else
    z=z_temp;
    w=(m(1,2)-m(2,1))/(4*z);
    x=(m(3,1)+m(1,3))/(4*z);
    y=(m(2,3)+m(3,2))/(4*z);
end
q=[w x y z];
    