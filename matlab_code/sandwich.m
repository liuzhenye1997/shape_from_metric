function result=sandwich(L,X,R)
XR=batch_quaternion_multiplication(X,R);
L(:,2:4)=-L(:,2:4);
result=batch_quaternion_multiplication(L,XR);