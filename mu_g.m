function mu_g=mu_g(info,P)
% P=P/6894.76; % Pa tp psi
mu_g = ((3e-10)*P.^2 + 1e-6.*P + 0.0133);%*10^-3; % cp tp Pa*s
end