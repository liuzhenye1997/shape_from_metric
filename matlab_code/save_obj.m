%�洢��name�ļ���
function  save_obj(vertex,face,name)
fid=fopen(name,'w');
vertex=vertex';
fprintf(fid,'v %.10f %.10f %.10f\r\n',vertex(:)); 
face=face';
fprintf(fid,'f %d %d %d\r\n',face(:)); 
fclose(fid);
end

