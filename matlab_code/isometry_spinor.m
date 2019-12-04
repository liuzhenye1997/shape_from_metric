%初始化face_psi,即λ
function lambda=isometry_spinor(face_number)
lambda=create_face_psi(face_number);
face_psi_norm=sqrt(sum(lambda.^2,2));
lambda=lambda./face_psi_norm;

%% 初始生成face_psi
function lambda=create_face_psi(face_number)
lambda(:,1)=rand(face_number,1);
lambda(:,2)=rand(face_number,1);
lambda(:,3)=rand(face_number,1);
lambda(:,4)=rand(face_number,1);
