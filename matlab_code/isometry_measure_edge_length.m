function length=isometry_measure_edge_length(vertex_src,vertex_dst,points)
length=sqrt(sum((points(vertex_dst',:)-points(vertex_src',:)).^2,2));