function Q_vector=qrotate(Q,vector)
vector=[zeros(size(vector,1),1) vector];
Q_transpose=[Q(:,1) -Q(:,2:4)];
Q_vector=batch_quaternion_multiplication(batch_quaternion_multiplication(Q,vector),Q_transpose);
Q_vector=Q_vector(:,2:4);