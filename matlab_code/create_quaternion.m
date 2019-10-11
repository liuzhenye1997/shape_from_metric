function Q=create_quaternion(theta,normal)
Q=[cos(theta/2) normal.*sin(theta/2)];