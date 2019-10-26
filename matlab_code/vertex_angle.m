%计算mesh的每一个面的三个角度
function angle=vertex_angle(length,vertex_prev,vertex_next)
l1=length';
l2=length(vertex_prev');
l2=l2(:)';
lopp=length(vertex_next');
lopp=lopp(:)';

c=(l1.^2+l2.^2-lopp.^2)./(2*l1.*l2);
angle=acos(c);