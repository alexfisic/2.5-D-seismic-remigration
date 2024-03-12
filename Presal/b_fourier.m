function [b]=b_fourier(n,T,f,x0,dx);

[ns nn]=size(f);
b=0;
omega=(2*pi)/T;
ni=round(x0/dx);
nf=round((x0+T)/dx);
nf=nf-1;

if(nf <= ns); 
  
for i=ni:nf;
    x=x0+(i-1)*dx;
    b=b+f(i,1)*sin(n*omega*x)*dx;  
end;

else;b=0;end;

b=b*(2/T);

