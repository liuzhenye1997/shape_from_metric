function quaternion=batch_angleaxis_to_quaternion(angleaxis)
norm=sqrt(sum(angleaxis.^2,2));
batch_number=size(angleaxis,1);
quaternion=zeros(batch_number,4);
angleaxis=angleaxis./norm;
quaternion(:,1)=1;
norm_index=norm>0;

quaternion(norm_index,:)=[cos(norm(norm_index)/2) angleaxis(norm_index,1).*sin(norm(norm_index)/2) angleaxis(norm_index,2).*sin(norm(norm_index)/2) angleaxis(norm_index,3).*sin(norm(norm_index)/2)];



    