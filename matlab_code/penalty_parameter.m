function [face_psi2im,face_psi2re,face_psi1re,face_psi1im]=penalty_parameter(face_psi1re_pre,face_psi2re_pre,face_psi2im_pre,face_psi1im_pre,face_psi2im,face_psi2re,face_psi1re,face_psi1im,a,dt)

if abs(a)<1*exp(-6)
    r=0;
else
    r=exp(-dt/a);
end

face_psi1im=face_psi1im*(1-r);
face_psi1re=face_psi1re*(1-r);
face_psi2re=face_psi2re*(1-r);
face_psi2im=face_psi2im*(1-r);

face_psi1im=face_psi1im+r*face_psi1im_pre;
face_psi2im=face_psi2im+r*face_psi2im_pre;
face_psi2re=face_psi2re+r*face_psi2re_pre;
face_psi1re=face_psi1re+r*face_psi1re_pre;

face_psi=[face_psi1im face_psi1re face_psi2im face_psi2re];
face_psi_norm=sqrt(sum(face_psi.^2,2));
face_psi1im=face_psi1im./face_psi_norm;
face_psi1re=face_psi1re./face_psi_norm;
face_psi2im=face_psi2im./face_psi_norm;
face_psi2re=face_psi2re./face_psi_norm;