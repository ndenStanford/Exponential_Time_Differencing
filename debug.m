clc
clear all
cd debug_plots
load('Implicit_883.4216_small_tol.mat')
max(max(delta_matrix_2D_S))
cd ..
% plotResultsImplicit
% figure
% p=[2900:1:4500];
% y=[];
% for i=1:length(p)
% y=[y,Bg(info,p(i))];
% end
% plot(p,y)
% xlabel('p')
% ylabel('Bg')
% % file=["30.mat","31.mat","32.mat"];
% file=["103.mat","104.mat","105.mat"];
% load('103.mat')
% Pressure_matrix_n_fix=Pressure_matrix_n;
% for i = 1:length(file)
% load(file(i))
% % Pressure_matrix_n(122/2)
% % an(122/2)
% % Pressure_matrix_n(61)
% x=exp_Lh*Pressure_matrix_n_fix;
% x(61)
% end

% figure
% % M2(M2<1e-4)=0;
% imagesc(M1)
% title('M1')
% colorbar