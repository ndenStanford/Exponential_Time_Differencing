function N_matrix_IMPES=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0)
[row,col]=size(Pressure_matrix_n);
P_matrix_n=zeros(2*row,col);
P_matrix_n(1:2:end-1,col)=Pressure_matrix_n;
P_matrix_n(2:2:end,col)=Saturation_matrix_n;
Pressure_matrix_n1=Pressure_matrix_n*(1+1e-5);
[D_matrix_out,d11_matrix,d12_matrix,d21_matrix,d22_matrix]=D_matrix_IMPES(info,P_matrix_n,Pressure_matrix_n1,delta_t_n);
Q_matrix_out=Q_matrix_IMPES(info,P_matrix_n,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix);

L=L_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n);
N_matrix_IMPES=(D_matrix_out*delta_t_n)\Q_matrix_out+(L-L0)*Pressure_matrix_n;
end
% spy(D_matrix)