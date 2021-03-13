% x = sort(unique(randi(50,1,10)));
% y = randi(10,1,length(x));
x=t_vec_0;
y=P_block;
xi = linspace(0,max(x),150);
yi = interp1(x,y,xi,'linear');
figure(1)
plot(x,y,'-*b')
hold on
plot(xi,yi,'pr')
hold off
grid
legend('Original','Interpolated')