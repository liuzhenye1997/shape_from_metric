%���ż�ߵĳ��ȡ���żͼ�ĵ�ȡ�������ġ�������ż�ߵĳ���Ϊ���ĵ���Ӧ�ߵľ��롣
function vertex_dual_length=dual_length_circumcentric(length,area,vertex_next,vertex_prev)
l0=length;
l1=length(vertex_next');
l2=length(vertex_prev');
l1=l1(:);
l2=l2(:);

area_temp=[area';area';area'];
area_temp=area_temp(:);
vertex_dual_length=(l1.*l1+l2.*l2-l0.*l0).*l0./(8*area_temp);