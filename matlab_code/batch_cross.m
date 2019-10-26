%p”Îq≤Ê≥À
function result=batch_cross(p,q)
result=[p(:,2).*q(:,3)-p(:,3).*q(:,2) p(:,3).*q(:,1)-p(:,1).*q(:,3) p(:,1).*q(:,2)-p(:,2).*q(:,1)];