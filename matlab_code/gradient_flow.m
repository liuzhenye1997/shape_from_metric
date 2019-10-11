function [points_psi1re,points_psi1im,points_psi2re,points_psi2im]=gradient_flow(L,dt,face_number,M,DMU,gradE_mu,points_psi1re,points_psi1im,points_psi2re,points_psi2im)
gradE_mu_temp=[gradE_mu(:,1);gradE_mu(:,2);gradE_mu(:,3);gradE_mu(:,4)];
gradE_2nd_term=(DMU')*(gradE_mu_temp(:));

area=M((size(M,1)+1)*(1:size(M,1))-size(M,1))';

Psi=[points_psi1re;points_psi1im;points_psi2re;points_psi2im];
Psi=Psi-dt.*(gradE_2nd_term./area);

LeftOp=M-dt.*L;
RightOp=M;
Psi=LeftOp\(RightOp*Psi);

points_psi1re=Psi(1:face_number);
points_psi1im=Psi(face_number+1:2*face_number);
points_psi2re=Psi(2*face_number+1:3*face_number);
points_psi2im=Psi(3*face_number+1:4*face_number);