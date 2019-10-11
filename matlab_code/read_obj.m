
function [vertex,face]=read_obj(name)
data=importdata(name);
vertex_number=0;
face_number=0;
for i=1:(size(data.textdata,1))  
    if char(data.textdata(i,1))=='v'
        vertex_number=vertex_number+1;
    elseif char(data.textdata(i,1))=='f'       
        face_number=face_number+1;
    end
end
vertex=zeros(vertex_number,3);
face=zeros(face_number,3);

sprintf('顶点个数为%d,面数为%d',vertex_number,face_number)  
vertex_no=1;
face_no=1;
size(vertex)
for i=1:(size(data.textdata,1))  
    if char(data.textdata(i,1))=='v'          
        vertex(vertex_no,:)=data.data(i,:);
        vertex_no=vertex_no+1;
    elseif char(data.textdata(i,1))=='f'
        face(face_no,:)=data.data(i,:);
        face_no=face_no+1;
    end
end   