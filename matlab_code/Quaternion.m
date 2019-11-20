classdef Quaternion
    %四元数    
    methods(Static = true)
        %利用旋转角度theta和单位旋转轴创建单位四元数      
        function  Q=create_quaternion(theta,normal)                           
            Q=[cos(theta/2) normal.*sin(theta/2)];          
        end
        
        %利用旋转矩阵创建四元数
        function Q=matrix_to_q(m)
            w_temp=sqrt(m(1,1)+m(2,2)+m(3,3)+1)/2;
            x_temp=sqrt(m(1,1)-m(2,2)-m(3,3)+1)/2;
            y_temp=sqrt(-m(1,1)+m(2,2)-m(3,3)+1)/2;
            z_temp=sqrt(-m(1,1)-m(2,2)+m(3,3)+1)/2;
            MAX=max([w_temp x_temp y_temp z_temp]);
            if w_temp==MAX
                w=w_temp;
                x=(m(2,3)-m(3,2))/(4*w);
                y=(m(3,1)-m(1,3))/(4*w);
                z=(m(1,2)-m(2,1))/(4*w);
            elseif x_temp==MAX
                x=x_temp;
                w=(m(2,3)-m(3,2))/(4*x);
                y=(m(1,2)+m(2,1))/(4*x);
                z=(m(3,1)+m(1,3))/(4*x);
            elseif y_temp==MAX
                y=y_temp;
                w=(m(3,1)-m(1,3))/(4*y);
                x=(m(1,2)+m(2,1))/(4*y);
                z=(m(2,3)+m(3,2))/(4*y);
            else
                z=z_temp;
                w=(m(1,2)-m(2,1))/(4*z);
                x=(m(3,1)+m(1,3))/(4*z);
                y=(m(2,3)+m(3,2))/(4*z);
            end
            Q=[w x y z];
        end       
        
        function Q=batch_matrix_to_q(m)
            w=sqrt(m(:,1)+m(:,5)+m(:,9)+1)/2;                    
            x=(m(:,6)-m(:,8))./(4*w);
            y=(m(:,7)-m(:,3))./(4*w);
            z=(m(:,2)-m(:,4))./(4*w);          
            Q=[w x y z];
        end         
        
        %只利用旋转轴创建四元数，其中旋转轴的模长为旋转的角度。输入的旋转轴只有一个
        function quaternion=angleaxis_to_q(angleaxis)
            norm=sqrt(sum(angleaxis.^2));
            if norm==0
                quaternion=[1 0 0 0];
            else
                angleaxis=angleaxis/norm;
                quaternion=[cos(norm/2) angleaxis(1)*sin(norm/2) angleaxis(2)*sin(norm/2) angleaxis(3)*sin(norm/2)];
            end
        end
        
        %只利用旋转轴创建四元数，其中旋转轴的模长为旋转的角度。输入的旋转轴的个数为第一维长度
        function quaternion=batch_angleaxis_to_q(angleaxis)
            norm=sqrt(sum(angleaxis.^2,2));
            batch_number=size(angleaxis,1);
            quaternion=zeros(batch_number,4);
            angleaxis=angleaxis./norm;
            quaternion(:,1)=1;
            norm_index=norm>0;
            quaternion(norm_index,:)=[cos(norm(norm_index)/2) angleaxis(norm_index,1).*sin(norm(norm_index)/2) angleaxis(norm_index,2).*sin(norm(norm_index)/2) angleaxis(norm_index,3).*sin(norm(norm_index)/2)];    
        end
        
        %两个四元数相乘
        function result=multiplication(p,q)
            result=[p(1)*q(1)-p(2)*q(2)-p(3)*q(3) p(1)*q(2)+p(2)*q(1)+p(3)*q(4)-p(4)*q(3) p(1)*q(3)-p(2)*q(4)+p(3)*q(1)+p(4)*q(2) p(1)*q(4)+p(2)*q(3)-p(3)*q(2)+p(4)*q(1)];    
        end
        
        %两批四元数相乘。输入的四元数的个数为第一维长度
        function result=batch_multiplication(p,q)
            result=[p(:,1).*q(:,1)-p(:,2).*q(:,2)-p(:,3).*q(:,3)-p(:,4).*q(:,4) ...
                p(:,1).*q(:,2)+p(:,2).*q(:,1)+p(:,3).*q(:,4)-p(:,4).*q(:,3) ...
                p(:,1).*q(:,3)-p(:,2).*q(:,4)+p(:,3).*q(:,1)+p(:,4).*q(:,2) ...
                p(:,1).*q(:,4)+p(:,2).*q(:,3)-p(:,3).*q(:,2)+p(:,4).*q(:,1)];
        end
        
        %向量vector在四元数的作用下旋转
        function Q_vector=qrotate(Q,vector)
            vector=[zeros(size(vector,1),1) vector];
            Q_transpose=[Q(:,1) -Q(:,2:4)];
            Q_vector=Quaternion.batch_multiplication(Quaternion.batch_multiplication(Q,vector),Q_transpose);
            Q_vector=Q_vector(:,2:4);
        end
        
        %单位四元数p的逆
        function result=qinvert(p)
            result=[p(:,1) -p(:,2) -p(:,3) -p(:,4)];
        end
        
        %计算从向量p旋转到向量q所需的四元数
        function quaternion=compute_dihedral(p,q)
            p_norm=sqrt(sum(p.^2,2));
            q_norm=sqrt(sum(q.^2,2));
            p=p./p_norm;
            q=q./q_norm;
            normal=[p(:,2).*q(:,3)-p(:,3).*q(:,2) p(:,3).*q(:,1)-p(:,1).*q(:,3) p(:,1).*q(:,2)-p(:,2).*q(:,1)];
            normal=normal./sqrt(sum(normal.^2,2));
            theta=acos(sum(p.*q,2));
            quaternion=[cos(theta/2) normal.*sin(theta/2)];
        end
        
    end
end

