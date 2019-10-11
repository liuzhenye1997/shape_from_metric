function vertex_theta=hedge_coordinate(angle,vertex_next,faces,face_hedge)
face_number=size(faces,1);
vertex_theta=zeros(size(faces))';
he0=face_hedge;
he=he0;
vertex_next_temp=vertex_next';
for i=1:face_number
    cumAngle=0;
    vertex_theta(he(i))=cumAngle;
    he_next=vertex_next_temp(he(i));
    cumAngle=cumAngle+angle(he_next)-pi;
    cumAngle = atan2(sin(cumAngle),cos(cumAngle));
    he(i)=he_next;

    while he(i)~=he0(i)
        vertex_theta(he(i))=cumAngle;
        he_next=vertex_next_temp(he(i));
        cumAngle=cumAngle+angle(he_next)-pi;
        cumAngle = atan2(sin(cumAngle),cos(cumAngle));
        he(i)=he_next;
    end
end
vertex_theta=vertex_theta';