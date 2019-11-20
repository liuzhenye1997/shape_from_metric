classdef Triangular_mesh
    %使用半边数据结构创建三角网格
    %这里point和vertex代指不同意义的点，即point是实际的点，而vertex是每个半边的源顶点。
    %故不同vertex在point里可能指向同一个点。在接下来的行文中用点代指point，顶点代指vertex。
    %以f开头的变量与面有关，与v开头的变量与顶点有关。以p开头的变量与点有关。
    properties
        angle          
        faces
        f_area
        f_number
        f_hedge%平面顶点所对应的半边
        f_ind%面的三条半边的序号
        length %每个面的三条边的长度
        points
        p_area
        p_number  
        v_face%顶点所在的面
        v_flip%对面的顶点
        v_dst %半边指向的点(point)
        v_next%半边指向的顶点
        v_prev%半边指向这个顶点的顶点
        v_src %半边的源点(point)
    end
    
    methods
        %创建时只有面和边的长度    
        function mesh = Triangular_mesh(faces,length)          
            mesh.f_area=face_area(length);         
            mesh.f_number = size(faces,1);
            mesh.f_ind=1:3*mesh.f_number;
            mesh.p_number=max(faces(:));
            mesh.points=zeros(size(length));%点的位置暂时用0代替            
            mesh.f_hedge=3*(0:mesh.f_number-1)+1;    
            mesh.faces=faces;
            length=length';
            mesh.length=length(:);
           
            %创建顶点属性
            mesh.v_face=[1:mesh.f_number;1:mesh.f_number;1:mesh.f_number]';
            mesh.v_src=mesh.faces;
            mesh.v_dst=[mesh.faces(:,2),mesh.faces(:,3),mesh.faces(:,1)];
            mesh.v_next=[(0:mesh.f_number-1)*3+2;(0:mesh.f_number-1)*3+3;(0:mesh.f_number-1)*3+1]';
            mesh.v_prev=[(0:mesh.f_number-1)*3+3;(0:mesh.f_number-1)*3+1;(0:mesh.f_number-1)*3+2]';           
            mesh.v_flip=vertex_flip(mesh.f_number,mesh.v_dst,mesh.v_src);          
            mesh.angle=vertex_angle(mesh.length,mesh.v_prev,mesh.v_next);  
        end
       
        function drawmesh(mesh)            
            drawmesh(mesh.faces, mesh.points);
        end   
        
    end
    methods(Static=true)
    end
end

