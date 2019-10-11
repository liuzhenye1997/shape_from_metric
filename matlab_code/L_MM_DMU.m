function [L,MM,DMU]=L_MM_DMU(k_r,k_b,k_c,k_t,e_r,e_b,e_c,e_t,face_number,star1,face_ind,face_src,face_dst,F,G,M,avgfac)
w11 = k_r.*e_r(:,1).*e_r(:,1) + k_b.*e_b(:,1).*e_b(:,1) + k_c.*e_c(:,1).*e_c(:,1) + k_t.*e_t(:,1).*e_t(:,1);
w12 = k_r.*e_r(:,1).*e_r(:,2) + k_b.*e_b(:,1).*e_b(:,2) + k_c.*e_c(:,1).*e_c(:,2) + k_t.*e_t(:,1).*e_t(:,2);
w13 = k_r.*e_r(:,1).*e_r(:,3) + k_b.*e_b(:,1).*e_b(:,3) + k_c.*e_c(:,1).*e_c(:,3) + k_t.*e_t(:,1).*e_t(:,3);
w14 = k_r.*e_r(:,1).*e_r(:,4) + k_b.*e_b(:,1).*e_b(:,4) + k_c.*e_c(:,1).*e_c(:,4) + k_t.*e_t(:,1).*e_t(:,4);
w22 = k_r.*e_r(:,2).*e_r(:,2) + k_b.*e_b(:,2).*e_b(:,2) + k_c.*e_c(:,2).*e_c(:,2) + k_t.*e_t(:,2).*e_t(:,2);
w23 = k_r.*e_r(:,2).*e_r(:,3) + k_b.*e_b(:,2).*e_b(:,3) + k_c.*e_c(:,2).*e_c(:,3) + k_t.*e_t(:,2).*e_t(:,3);
w24 = k_r.*e_r(:,2).*e_r(:,4) + k_b.*e_b(:,2).*e_b(:,4) + k_c.*e_c(:,2).*e_c(:,4) + k_t.*e_t(:,2).*e_t(:,4);
w33 = k_r.*e_r(:,3).*e_r(:,3) + k_b.*e_b(:,3).*e_b(:,3) + k_c.*e_c(:,3).*e_c(:,3) + k_t.*e_t(:,3).*e_t(:,3);
w34 = k_r.*e_r(:,3).*e_r(:,4) + k_b.*e_b(:,3).*e_b(:,4) + k_c.*e_c(:,3).*e_c(:,4) + k_t.*e_t(:,3).*e_t(:,4);
w44 = k_r.*e_r(:,4).*e_r(:,4) + k_b.*e_b(:,4).*e_b(:,4) + k_c.*e_c(:,4).*e_c(:,4) + k_t.*e_t(:,4).*e_t(:,4);

W11=sparse(1:3*face_number,1:3*face_number,w11)*star1;
W12=sparse(1:3*face_number,1:3*face_number,w12)*star1;
W13=sparse(1:3*face_number,1:3*face_number,w13)*star1;
W14=sparse(1:3*face_number,1:3*face_number,w14)*star1;
W22=sparse(1:3*face_number,1:3*face_number,w22)*star1;
W23=sparse(1:3*face_number,1:3*face_number,w23)*star1;
W24=sparse(1:3*face_number,1:3*face_number,w24)*star1;
W33=sparse(1:3*face_number,1:3*face_number,w33)*star1;
W34=sparse(1:3*face_number,1:3*face_number,w34)*star1;
W44=sparse(1:3*face_number,1:3*face_number,w44)*star1;

W1=[W11,W12,W13,W14];
W2=[W12,W22,W23,W24];
W3=[W13,W23,W33,W34];
W4=[W14,W24,W34,W44];
W=[W1;W2;W3;W4];

face_src=face_src';
face_dst=face_dst';
d_row=[face_ind,face_ind]';
d_col=[face_src(:);face_dst(:)];

D1_val=[-ones(3*face_number,1),G(:,1)];
D2_val=[zeros(3*face_number,1),G(:,2)];
D3_val=[zeros(3*face_number,1),G(:,3)];
D4_val=[zeros(3*face_number,1),G(:,4)];


D1=sparse(d_row,d_col,D1_val,3*face_number,face_number);
D2=sparse(d_row,d_col,D2_val,3*face_number,face_number);
D3=sparse(d_row,d_col,D3_val,3*face_number,face_number);
D4=sparse(d_row,d_col,D4_val,3*face_number,face_number);

D_1=[D1,-D2,-D3,-D4];
D_2=[D2,D1,-D4,D3];
D_3=[D3,D4,D1,-D2];
D_4=[D4,-D3,D2,D1];
D=[D_1;D_2;D_3;D_4];


L=-D'*W*D;
ZM=sparse(size(M,1),size(M,2));
MN_1=[M,ZM,ZM,ZM];
MN_2=[ZM,M,ZM,ZM];
MN_3=[ZM,ZM,M,ZM];
MN_4=[ZM,ZM,ZM,M];
MM=[MN_1;MN_2;MN_3;MN_4];

AVG1_val=[avgfac;(F(:,1).*avgfac)];
AVG2_val=[zeros(3*face_number,1);(F(:,2).*avgfac)];
AVG3_val=[zeros(3*face_number,1);(F(:,3).*avgfac)];
AVG4_val=[zeros(3*face_number,1);(F(:,4).*avgfac)];

AVG1=sparse(d_row,d_col,AVG1_val,3*face_number,face_number);
AVG2=sparse(d_row,d_col,AVG2_val,3*face_number,face_number);
AVG3=sparse(d_row,d_col,AVG3_val,3*face_number,face_number);
AVG4=sparse(d_row,d_col,AVG4_val,3*face_number,face_number);
AVG_1=[AVG1,-AVG2,-AVG3,-AVG4];
AVG_2=[AVG2,AVG1,-AVG4,AVG3];
AVG_3=[AVG3,AVG4,AVG1,-AVG2];
AVG_4=[AVG4,-AVG3,AVG2,AVG1];
AVG=[AVG_1;AVG_2;AVG_3;AVG_4];

P11=sparse(1:3*face_number,1:3*face_number,1-e_r(:,1).*e_r(:,1));
P12=sparse(1:3*face_number,1:3*face_number,-e_r(:,1).*e_r(:,2));
P13=sparse(1:3*face_number,1:3*face_number,-e_r(:,1).*e_r(:,3));
P14=sparse(1:3*face_number,1:3*face_number,-e_r(:,1).*e_r(:,4));
P22=sparse(1:3*face_number,1:3*face_number,1-e_r(:,2).*e_r(:,2));
P23=sparse(1:3*face_number,1:3*face_number,-e_r(:,2).*e_r(:,3));
P24=sparse(1:3*face_number,1:3*face_number,-e_r(:,2).*e_r(:,4));
P33=sparse(1:3*face_number,1:3*face_number,1-e_r(:,3).*e_r(:,3));
P34=sparse(1:3*face_number,1:3*face_number,-e_r(:,3).*e_r(:,4));
P44=sparse(1:3*face_number,1:3*face_number,1-e_r(:,4).*e_r(:,4));
P_1=[P11,P12,P13,P14];
P_2=[P12,P22,P23,P24];
P_3=[P13,P23,P33,P34];
P_4=[P14,P24,P34,P44];
P=[P_1;P_2;P_3;P_4];

DMU=P*AVG;


