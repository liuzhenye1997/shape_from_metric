function [Prev_k_real,Prev_k_connection,Prev_k_torsion,Prev_k_bending]=create_k(iteration,face_number)
if iteration<10
    Prev_k_real=ones(3*face_number,1);
    Prev_k_connection=ones(3*face_number,1);
    Prev_k_torsion=ones(3*face_number,1);
    Prev_k_bending=ones(3*face_number,1);
elseif iteration<20
    Prev_k_real=ones(3*face_number,1);
    Prev_k_connection=ones(3*face_number,1);
    Prev_k_torsion=ones(3*face_number,1)*0.1;
    Prev_k_bending=ones(3*face_number,1)*0.1;
elseif iteration<50
    Prev_k_real=ones(3*face_number,1);
    Prev_k_connection=ones(3*face_number,1);
    Prev_k_torsion=ones(3*face_number,1)*0.1;
    Prev_k_bending=ones(3*face_number,1)*0.01;
else
    Prev_k_real=ones(3*face_number,1);
    Prev_k_connection=ones(3*face_number,1);
    Prev_k_torsion=ones(3*face_number,1)*0.1;
    Prev_k_bending=zeros(3*face_number,1);
end