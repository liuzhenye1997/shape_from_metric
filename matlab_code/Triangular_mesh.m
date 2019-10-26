classdef Triangular_mesh
    %ʹ�ð�����ݽṹ������������
    %����point��vertex��ָ��ͬ����ĵ㣬��point��ʵ�ʵĵ㣬��vertex��ÿ����ߵ�Դ���㡣
    %�ʲ�ͬvertex��point�����ָ��ͬһ���㡣�ڽ��������������õ��ָpoint�������ָvertex��
    %��f��ͷ�ı��������йأ���v��ͷ�ı����붥���йء���p��ͷ�ı�������йء�
    properties
        angle          
        faces
        f_area
        f_number
        f_hedge%ƽ�涥������Ӧ�İ��
        f_ind%���������ߵ����
        length %ÿ����������ߵĳ���
        points
        p_area
        p_number  
        v_face%�������ڵ���
        v_flip%����Ķ���
        v_dst %���ָ��ĵ�(point)
        v_next%���ָ��Ķ���
        v_prev%���ָ���������Ķ���
        v_src %��ߵ�Դ��(point)
    end
    
    methods
        %����ʱֻ����ͱߵĳ���    
        function mesh = Triangular_mesh(faces,length)
            
            mesh.f_area=face_area(length);         
            mesh.f_number = size(faces,1);
            mesh.f_ind=1:3*mesh.f_number;
            mesh.p_number=max(faces(:));
            mesh.points=zeros(size(length));%���λ����ʱ��0����            
            mesh.f_hedge=3*(0:mesh.f_number-1)+1;    
            mesh.faces=faces;
            length=length';
            mesh.length=length(:);
            
            %������������
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

