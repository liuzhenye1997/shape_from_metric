function gradE_mu=grad_energy_wrt_mu(xi,Z,I,IZ,mu,face_weight)
hp_1_xi=xi;
hp_Z_xi=batch_quaternion_multiplication(qinvert(Z),xi);
hp_I_xi=batch_quaternion_multiplication(qinvert(I),xi);
hp_IZ_xi=batch_quaternion_multiplication(qinvert(IZ),xi);

ip_mu_xi=batch_quaternion_multiplication(qinvert(mu),hp_1_xi);
ip_Zmu_xi=batch_quaternion_multiplication(qinvert(mu),hp_Z_xi);
ip_Imu_xi=batch_quaternion_multiplication(qinvert(mu),hp_I_xi);
ip_IZmu_xi=batch_quaternion_multiplication(qinvert(mu),hp_IZ_xi);


gradE_mu=face_weight.*(hp_1_xi.*ip_mu_xi(:,1)+hp_Z_xi.*ip_Zmu_xi(:,1)+hp_I_xi.*ip_Imu_xi(:,1)+hp_IZ_xi.*ip_IZmu_xi(:,1));
