%���ñ߳������
function area=face_area(length)
p=sum(length,2)/2;
area=sqrt(p.*(p-length(:,1)).*(p-length(:,2)).*(p-length(:,3)));

