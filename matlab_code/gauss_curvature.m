function points_gauss_curvature=gauss_curvature(faces,angle,point_number)
points_gauss_curvature=zeros(1,point_number);
face_number=size(faces,1);
for i=1:face_number
    points_gauss_curvature(faces(i,1))=points_gauss_curvature(faces(i,1))+angle(3*i-2);
    points_gauss_curvature(faces(i,2))=points_gauss_curvature(faces(i,2))+angle(3*i);
    points_gauss_curvature(faces(i,3))=points_gauss_curvature(faces(i,3))+angle(3*i-1);
end