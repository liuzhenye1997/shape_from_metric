function points_isBad=label_bad_vertices(points_holonomyFromA,points_holonomy,type)
point_number=size(points_holonomy,2);
points_isBad=zeros(point_number,1);
if type==1
    for i=1:point_number
        evenOrOdd=floor((points_holonomy(i)-points_holonomyFromA(i))/pi+0.5);
        points_isBad(i)=mod(evenOrOdd,2);
    end
else
    for i=1:point_number
        plusorminus=cos(points_holonomyFromA(i))*cos(points_holonomy(i))+sin(points_holonomyFromA(i))*sin(points_holonomy(i));
        points_isBad(i)=plusorminus<0;
    end
end
