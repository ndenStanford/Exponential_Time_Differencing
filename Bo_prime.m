function Bo_prime=Bo_prime(info,P)
% % P=P/6894.76; % Pa to psi 
% P_bub=3400; % psi 
% P_atm=14.7; % psi 
% co=0.8e-5;
if P<info.P_bub
Bo_prime=(8e-5)*exp(-8e-5*(info.P_atm-P));  
else     
Bo_prime=(-info.co)*exp(-8e-5*(info.P_atm-info.P_bub))*exp(-info.co*(P-info.P_bub));  
end
end