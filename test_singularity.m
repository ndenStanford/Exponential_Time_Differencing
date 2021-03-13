n=9;
A=diag(2*ones(n,1))+diag(-1*ones(n-1,1),1)+diag(-1*ones(n-1,1),-1);
A(1,1)=1;
A(n,n)=1;