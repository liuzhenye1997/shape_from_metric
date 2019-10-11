function quaternion=angleaxis_to_quaternion(angleaxis)
norm=sqrt(sum(angleaxis.^2));
if norm==0
    quaternion=[1 0 0 0];
else
    angleaxis=angleaxis/norm;
    quaternion=[cos(norm/2) angleaxis(1)*sin(norm/2) angleaxis(2)*sin(norm/2) angleaxis(3)*sin(norm/2)];
end