%³õÊ¼Éú³Éface_psi
function [face_psi1re,face_psi1im,face_psi2re,face_psi2im]=create_face_psi(seed,face_psi1re,face_psi1im,face_psi2re,face_psi2im)
face_number=size(face_psi1re,1);
face_psi1re(:)=sin(1:face_number+seed);
face_psi1im(:)=sin((1:face_number)*pi*pi+seed);
face_psi2re(:)=sin((1:face_number)*pi*exp(1)+seed);
face_psi2im(:)=sin((1:face_number)*exp(1)+seed);
