function Rs_prime=Rs_prime(info,P)
% P_bub=3400; % psi 
% % P=P/6894.76; % Pa to psi
if P>info.P_bub
Rs_prime=0;    
else    
Rs_prime=(178.11^2/5.615)*1.3*P^0.3*(1/info.P_bub)^1.3;
end
end