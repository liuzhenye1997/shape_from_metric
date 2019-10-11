function root_candidate=find_root_candidate(points_local_isBdy)
point_number=size(points_local_isBdy,1);
root_candidate=1;
for i=1:point_number
    isBdy=points_local_isBdy(i);
    if isBdy==1
        root_candidate=i;
        break;
    end
end
