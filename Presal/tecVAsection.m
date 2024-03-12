function secVAout=tecVAsection(nt,nx,U);

U2=zeros(nt,nx);
U3=zeros(nt,nx);

tec=zeros(nt,nx);

M=16;

ii=complex(0,1);

for ix=1:nx;

for it=1:nt;
    inf = it-(M/2); sup = it+(M/2);
    if(inf > 0 && sup <= nt);
    for i=inf:sup;
    tec(it,ix) = tec(it,ix) + U(i,ix)* U(i,ix);
    end;
    tec(it,ix)=sqrt(tec(it,ix)/M);
    else; end;
end;

end;

for ix=1:nx;
    
    U2(:,ix)=fft(tec(:,ix));
    
    for it = 1:round(nt/2 + 2);
    U2(it,ix) = -ii*U2(it,ix);
    end;

    for it = round(nt/2 + 3):nt;
    U2(it,ix) = ii*U2(it,ix);
    end;

end;

for ix=1:nx;
    U3(:,ix)=real(ifft(U2(:,ix)));
end;

%for i=1:nx;
%U3(:,i)=U3(:,i)/max(abs(U3(:,i)));
%end;

secVAout=U3;

clear U2 U3 tec;



