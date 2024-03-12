function [taper]=tuckey(nt,alpha);

for n=2:nt-1;
    if(0 <= n & n <= (alpha*nt/2));
        taper(1,n) = 0.5*(1 + cos(pi*((2*n/(alpha*nt))- 1)));
    elseif((alpha*nt/2) <= n & n <= nt*(1-alpha/2));
        taper(1,n) = 1.0;
    elseif (nt*(1-alpha/2) <= n & n <= nt);
        taper(1,n) = 0.5*(1 + cos(pi*((2*n/(alpha*nt))- (2/alpha) + 1)));
    else; end;
end;

taper(1,1)=0.0;
taper(1,nt)=0.0;
