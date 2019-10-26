classdef Quaternion
    %��Ԫ��    
    methods(Static = true)
        %������ת�Ƕ�theta�͵�λ��ת�ᴴ����λ��Ԫ��      
        function  Q=create_quaternion(theta,normal)                           
            Q=[cos(theta/2) normal.*sin(theta/2)];          
        end
        %������ת���󴴽���Ԫ��
        function Q=matrix_to_quaternion(m)
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
        %ֻ������ת�ᴴ����Ԫ����������ת���ģ��Ϊ��ת�ĽǶȡ��������ת��ֻ��һ��
        function quaternion=angleaxis_to_quaternion(angleaxis)
            norm=sqrt(sum(angleaxis.^2));
            if norm==0
                quaternion=[1 0 0 0];
            else
                angleaxis=angleaxis/norm;
                quaternion=[cos(norm/2) angleaxis(1)*sin(norm/2) angleaxis(2)*sin(norm/2) angleaxis(3)*sin(norm/2)];
            end
        end
        %ֻ������ת�ᴴ����Ԫ����������ת���ģ��Ϊ��ת�ĽǶȡ��������ת��ĸ���Ϊ��һά����
        function quaternion=batch_angleaxis_to_quaternion(angleaxis)
            norm=sqrt(sum(angleaxis.^2,2));
            batch_number=size(angleaxis,1);
            quaternion=zeros(batch_number,4);
            angleaxis=angleaxis./norm;
            quaternion(:,1)=1;
            norm_index=norm>0;

            quaternion(norm_index,:)=[cos(norm(norm_index)/2) angleaxis(norm_index,1).*sin(norm(norm_index)/2) angleaxis(norm_index,2).*sin(norm(norm_index)/2) angleaxis(norm_index,3).*sin(norm(norm_index)/2)];    
        end
        %������Ԫ�����
        function result=quaternion_multiplication(p,q)
            result=[p(1)*q(1)-p(2)*q(2)-p(3)*q(3) p(1)*q(2)+p(2)*q(1)+p(3)*q(4)-p(4)*q(3) p(1)*q(3)-p(2)*q(4)+p(3)*q(1)+p(4)*q(2) p(1)*q(4)+p(2)*q(3)-p(3)*q(2)+p(4)*q(1)];    
        end
        %������Ԫ����ˡ��������Ԫ���ĸ���Ϊ��һά����
        function result=batch_quaternion_multiplication(p,q)
            result=[p(:,1).*q(:,1)-p(:,2).*q(:,2)-p(:,3).*q(:,3)-p(:,4).*q(:,4) ...
                p(:,1).*q(:,2)+p(:,2).*q(:,1)+p(:,3).*q(:,4)-p(:,4).*q(:,3) ...
                p(:,1).*q(:,3)-p(:,2).*q(:,4)+p(:,3).*q(:,1)+p(:,4).*q(:,2) ...
                p(:,1).*q(:,4)+p(:,2).*q(:,3)-p(:,3).*q(:,2)+p(:,4).*q(:,1)];
        end
        %����vector����Ԫ������������ת
        function Q_vector=qrotate(Q,vector)
            vector=[zeros(size(vector,1),1) vector];
            Q_transpose=[Q(:,1) -Q(:,2:4)];
            Q_vector=Quaternion.batch_quaternion_multiplication(Quaternion.batch_quaternion_multiplication(Q,vector),Q_transpose);
            Q_vector=Q_vector(:,2:4);
        end
        %��λ��Ԫ��p����
        function result=qinvert(p)
            result=[p(:,1) -p(:,2) -p(:,3) -p(:,4)];
        end
        %���������p��ת������q�������Ԫ��
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

