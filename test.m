clear all
close all
clc
load('plot_ETDRK4_low_0.01.mat')
P_block_low=P_block;
load('plot_RK4_mid_0.01.mat')
P_block_mid=P_block;
diff=abs(P_block_low-P_block_mid);
% there is no etdrk4-Lie mid and loe