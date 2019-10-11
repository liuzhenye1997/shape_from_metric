function vertex_A=define_parallel_transport_angle_A(vertex_theta,vertex_flip)
vertex_theta_temp=vertex_theta';
vertex_A=zeros(size(vertex_theta));
face_number=size(vertex_theta,1);
for i=1:face_number
    for j=1:3
        if vertex_flip(3*(i-1)+j)~=-1
            thetaOpp=vertex_theta_temp(vertex_flip(3*(i-1)+j));
            vertex_A(i,j)=thetaOpp-vertex_theta(i,j)-pi;
        end
        vertex_A(i,j)=atan2(sin(vertex_A(i,j)),cos(vertex_A(i,j)));
    end
end