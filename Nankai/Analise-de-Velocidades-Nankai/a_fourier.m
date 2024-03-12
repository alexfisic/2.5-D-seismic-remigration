function [a,ni,nf]=a_fourier(n,T,f,x0,dx);

[ns nn]=size(f); 
a=0;
omega=(2*pi)/T;
ni=round(x0/dx);
nf=round((x0+(1*T))/dx);
nf=nf-1;

if(nf <= ns);

for i=ni:nf;
    x=x0+(i-1)*dx;
    a=a+f(i,1)*cos(n*omega*x)*dx;   
end;

else;a=0;end;

a=a*(2/T);




