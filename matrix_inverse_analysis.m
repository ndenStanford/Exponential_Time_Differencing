clc
clear all
close all
load('L_matrix_IMPES_vars')
tic
T_matrix_out_inv_1=inv(T_matrix_out);
test_1=T_matrix_out_inv_1*T_matrix_out;
toc

tic
T_matrix_out_inv_2=T_matrix_out\eye(529);
test_2=T_matrix_out_inv_2*T_matrix_out;
toc

T_matrix_out_inv_diff=T_matrix_out_inv_1-T_matrix_out_inv_2;
y=rand(529,1);
transpose(y)*T_matrix_out*y
eigs_val=eigs(T_matrix_out);