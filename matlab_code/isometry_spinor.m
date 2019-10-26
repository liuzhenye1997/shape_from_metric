%≥ı ºªØface_psi,º¥¶À
function [face_psi1re,face_psi1im,face_psi2re,face_psi2im]=isometry_spinor(face_number)
face_psi1re=zeros(face_number,1);
face_psi1im=zeros(face_number,1);
face_psi2re=zeros(face_number,1);
face_psi2im=zeros(face_number,1);

[face_psi1re,face_psi1im,face_psi2re,face_psi2im]=create_face_psi(0,face_psi1re,face_psi1im,face_psi2re,face_psi2im);
% face_psi1re(:)=1:face_number;
% face_psi1im(:)=2*(2:face_number++1);
% face_psi2re(:)=3*(3:face_number+2);
% face_psi2im(:)=4*(4:face_number+3);
face_psi=[face_psi1re face_psi1im face_psi2re face_psi2im];
face_psi_norm=sqrt(sum(face_psi.^2,2));
face_psi1re=face_psi1re./face_psi_norm;
face_psi1im=face_psi1im./face_psi_norm;
face_psi2re=face_psi2re./face_psi_norm;
face_psi2im=face_psi2im./face_psi_norm;