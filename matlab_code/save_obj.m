%将结果输入new_Bunny_head.obj
function  save_obj(vertex,face,name)
fid=fopen([name],'w');
x=1;
for i=1:size(vertex,1)
        fprintf(fid,'v %.10f %.10f %.10f\r\n',vertex(i,1),vertex(i,2),vertex(i,3)); 
end
for i=1:size(face,1)
    fprintf(fid,'f %d %d %d\r\n',face(i,1),face(i,2),face(i,3));      
end
fclose(fid);
end