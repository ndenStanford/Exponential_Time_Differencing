function L_matrix_IMPES=L_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n)
[row,col]=size(Pressure_matrix_n);
P_matrix_n=zeros(2*row,col);
P_matrix_n(1:2:end-1)=Pressure_matrix_n;
P_matrix_n(2:2:end)=Saturation_matrix_n;
Pressure_matrix_n1=Pressure_matrix_n*(1+1e-5);

[D_matrix_out,d11_matrix,d12_matrix,d21_matrix,d22_matrix]=D_matrix_IMPES(info,P_matrix_n,Pressure_matrix_n1,delta_t_n);
T_matrix_out=T_matrix_IMPES(info,P_matrix_n,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
L_matrix_IMPES=(D_matrix_out*delta_t_n)\T_matrix_out;
%L_matrix_IMPES_inv=T_matrix_out\(D_matrix_out*delta_t_n);
end