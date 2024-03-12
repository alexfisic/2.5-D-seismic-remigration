close all;
index=input('Posicao de referencia (0-800)-->');
load urucu_velo_3.dat;
velo_ml=urucu_velo_3';
v=velo_ml(:,index)*1000;
v=v';
z=[10:10:3500];
%v=[1750 1750 2300 3750 5900 3850 6150 4000 5300 5500 6150 5450 6150 5600 5400 5200 5600 5400 5200 6150 4850 4650 6100 5250];
%z=[0:152.1739:3500];
dz=152.1739;
v2=1750+.8593*z;

X=z;
Y=v;
XX=X.*X;
XXX=X.*XX;
XXXX=X.*XXX;
XY=X.*Y;
XXY=XX.*Y;
YY=Y.*Y;
X5=XXXX.*X;
X6=X5.*X;
X7=X6.*X;
X8=X7.*X;
X9=X8.*X;
X10=X9.*X;
N=length(v);

a0=(sum(Y)*sum(XX)-sum(X)*sum(XY))/(length(v)*sum(XX)-sum(X)*sum(X));
a1=(length(v)*sum(XY)-sum(X)*sum(Y))/(length(v)*sum(XX)-sum(X)*sum(X));

v3=a0+a1*z;

figure;plot(z,v,'gV',z,v,'-.',z,v2,z,v3,'c--');grid;


v_medio=zeros(size(v));


for i=1:length(z);
    z1=z(1)+(i-1)*dz;
    nz=round(z1/dz);
    for j=1:nz;
    v_medio(i)=v_medio(i)+z1*v(j)*dz;
    end;
end;

v_medio=v_medio/length(z);

%figure;plot(z,v,'gV',z,v_medio);

P=[N sum(X) sum(XX) sum(XXX) sum(XXXX) sum(X5); sum(X) sum(XX) sum(XXX) sum(XXXX) sum(X5) sum(X6); sum(XX) sum(XXX) sum(XXXX) sum(XXXX.*X) sum(X6) sum(X7); sum(XXX) sum(XXXX) sum(XXXX.*X) sum(XXXX.*XX) sum(X7) sum(X8); sum(XXXX) sum(X5) sum(X6) sum(X7) sum(X8) sum(X9); sum(X5) sum(X6) sum(X7) sum(X8) sum(X9) sum(X10)];
det_a=[sum(Y) sum(X) sum(XX) sum(XXX) sum(XXXX) sum(X5); sum(XY) sum(XX) sum(XXX) sum(XXXX) sum(X5) sum(X6); sum(XXY) sum(XXX) sum(XXXX) sum(XXXX.*X) sum(X6) sum(X7); sum(XXX.*Y) sum(XXXX) sum(XXXX.*X) sum(XXXX.*XX) sum(X7) sum(X8);sum(XXXX.*Y) sum(X5) sum(X6) sum(X7) sum(X8) sum(X9);sum(X5.*Y) sum(X6) sum(X7) sum(X8) sum(X9) sum(X10)];
det_b=[N sum(Y) sum(XX) sum(XXX) sum(XXXX) sum(X5); sum(X) sum(XY) sum(XXX) sum(XXXX) sum(X5) sum(X6); sum(XX) sum(XXY) sum(XXXX) sum(XXXX.*X) sum(X6) sum(X7); sum(XXX) sum(XXX.*Y) sum(XXXX.*X) sum(XXXX.*XX) sum(X7) sum(X8); sum(XXXX) sum(XXXX.*Y) sum(X6) sum(X7) sum(X8) sum(X9); sum(X5) sum(X5.*Y) sum(X7) sum(X8) sum(X9) sum(X10)];
det_c=[N sum(X) sum(Y) sum(XXX) sum(XXXX) sum(X5); sum(X) sum(XX) sum(XY) sum(XXXX) sum(X5) sum(X6); sum(XX) sum(XXX) sum(XXY) sum(XXXX.*X) sum(X6) sum(X7); sum(XXX) sum(XXXX) sum(XXX.*Y) sum(XXXX.*XX) sum(X7) sum(X8); sum(XXXX) sum(X5) sum(XXXX.*Y) sum(X7) sum(X8) sum(X9); sum(X5) sum(X6) sum(X5.*Y) sum(X8) sum(X9) sum(X10)];
det_d=[N sum(X) sum(XX) sum(Y) sum(XXXX) sum(X5); sum(X) sum(XX) sum(XXX) sum(XY) sum(X5) sum(X6); sum(XX) sum(XXX) sum(XXXX) sum(XXY) sum(X6) sum(X7); sum(XXX) sum(XXXX) sum(XXXX.*X) sum(XXX.*Y) sum(X7) sum(X8); sum(XXXX) sum(X5) sum(X6) sum(XXXX.*Y) sum(X8) sum(X9); sum(X5) sum(X6) sum(X7) sum(X5.*Y) sum(X9)  sum(X10)];
det_e=[N sum(X) sum(XX) sum(XXX) sum(Y) sum(X5); sum(X) sum(XX) sum(XXX) sum(XXXX) sum(XY) sum(X6); sum(XX) sum(XXX) sum(XXXX) sum(X5) sum(XXY) sum(X7); sum(XXX) sum(XXXX) sum(X5) sum(X6) sum(XXX.*Y) sum(X8); sum(XXXX) sum(X5) sum(X6) sum(X7) sum(XXXX.*Y) sum(X9);sum(X5) sum(X6) sum(X7) sum(X8) sum(X5.*Y) sum(X10)];
det_f=[N sum(X) sum(XX) sum(XXX) sum(XXXX) sum(Y); sum(X) sum(XX) sum(XXX) sum(XXXX) sum(X5) sum(XY); sum(XX) sum(XXX) sum(XXXX) sum(X5) sum(X6) sum(XXY); sum(XXX) sum(XXXX) sum(X5) sum(X6) sum(X7) sum(XXX.*Y); sum(XXXX) sum(X5) sum(X6) sum(X7) sum(X8) sum(XXXX.*Y); sum(X5) sum(X6) sum(X7) sum(X8) sum(X9) sum(X5.*Y)];

a3=det(det_a)/det(P);
a4=det(det_b)/det(P);
a5=det(det_c)/det(P);
a6=det(det_d)/det(P);
a7=det(det_e)/det(P);
a8=det(det_f)/det(P);

v4=a3+a4.*z+a5.*z.*z;
v5=a3+a4.*z+a5.*z.*z+a6.*z.*z.*z;
v6=a3+a4.*z+a5.*z.*z+a6.*z.*z.*z+a7.*z.*z.*z.*z;
v7=a3+a4.*z+a5.*z.*z+a6.*z.*z.*z+a7.*z.*z.*z.*z+a8.*z.*z.*z.*z.*z;

figure, plot(z,v,'gV',z,v,'-.',z,v5,'kO',z,v6,'r+',z,v7,'b*');grid;

figure,plot(z,abs(((v-v6)./v)*100),'r-.');grid;

omega=(2*pi/3500);

c0=a_fourier(0,3500,v,10,10);
c1=a_fourier(1,3500,v,10,10);
c2=a_fourier(2,3500,v,10,10);
c3=a_fourier(3,3500,v,10,10);
c4=a_fourier(4,3500,v,10,10);
c5=a_fourier(5,3500,v,10,10);
c6=a_fourier(6,3500,v,10,10);
c7=a_fourier(7,3500,v,10,10);
c8=a_fourier(8,3500,v,10,10);
c9=a_fourier(9,3500,v,10,10);
c10=a_fourier(10,3500,v,10,10);
c11=a_fourier(11,3500,v,10,10);
c12=a_fourier(12,3500,v,10,10);
c13=a_fourier(13,3500,v,10,10);
c14=a_fourier(14,3500,v,10,10);
c15=a_fourier(15,3500,v,10,10);
c16=a_fourier(16,3500,v,10,10);
c17=a_fourier(17,3500,v,10,10);
c18=a_fourier(18,3500,v,10,10);
c19=a_fourier(19,3500,v,10,10);
c20=a_fourier(20,3500,v,10,10);
c21=a_fourier(21,3500,v,10,10);
c22=a_fourier(22,3500,v,10,10);
c23=a_fourier(23,3500,v,10,10);
c24=a_fourier(24,3500,v,10,10);
c25=a_fourier(25,3500,v,10,10);
c26=a_fourier(26,3500,v,10,10);
c27=a_fourier(27,3500,v,10,10);
c28=a_fourier(28,3500,v,10,10);
c29=a_fourier(29,3500,v,10,10);
c30=a_fourier(30,3500,v,10,10);
c31=a_fourier(31,3500,v,10,10);
c32=a_fourier(32,3500,v,10,10);
c33=a_fourier(33,3500,v,10,10);
c34=a_fourier(34,3500,v,10,10);
c35=a_fourier(35,3500,v,10,10);
c36=a_fourier(36,3500,v,10,10);
c37=a_fourier(37,3500,v,10,10);
c38=a_fourier(38,3500,v,10,10);
c39=a_fourier(39,3500,v,10,10);
c40=a_fourier(40,3500,v,10,10);


d1=b_fourier(1,3500,v,10,10);
d2=b_fourier(2,3500,v,10,10);
d3=b_fourier(3,3500,v,10,10);
d4=b_fourier(4,3500,v,10,10);
d5=b_fourier(5,3500,v,10,10);
d6=b_fourier(6,3500,v,10,10);
d7=b_fourier(7,3500,v,10,10);
d8=b_fourier(8,3500,v,10,10);
d9=b_fourier(9,3500,v,10,10);
d10=b_fourier(10,3500,v,10,10);
d11=b_fourier(11,3500,v,10,10);
d12=b_fourier(12,3500,v,10,10);
d13=b_fourier(13,3500,v,10,10);
d14=b_fourier(14,3500,v,10,10);
d15=b_fourier(15,3500,v,10,10);
d16=b_fourier(16,3500,v,10,10);
d17=b_fourier(17,3500,v,10,10);
d18=b_fourier(18,3500,v,10,10);
d19=b_fourier(19,3500,v,10,10);
d20=b_fourier(20,3500,v,10,10);
d21=b_fourier(21,3500,v,10,10);
d22=b_fourier(22,3500,v,10,10);
d23=b_fourier(23,3500,v,10,10);
d24=b_fourier(24,3500,v,10,10);
d25=b_fourier(25,3500,v,10,10);
d26=b_fourier(26,3500,v,10,10);
d27=b_fourier(27,3500,v,10,10);
d28=b_fourier(28,3500,v,10,10);
d29=b_fourier(29,3500,v,10,10);
d30=b_fourier(30,3500,v,10,10);
d31=b_fourier(31,3500,v,10,10);
d32=b_fourier(32,3500,v,10,10);
d33=b_fourier(33,3500,v,10,10);
d34=b_fourier(34,3500,v,10,10);
d35=b_fourier(35,3500,v,10,10);
d36=b_fourier(36,3500,v,10,10);
d37=b_fourier(37,3500,v,10,10);
d38=b_fourier(38,3500,v,10,10);
d39=b_fourier(39,3500,v,10,10);
d40=b_fourier(40,3500,v,10,10);


v7=c0+c1*cos(omega*z)+c2*cos(2*omega*z)+c3*cos(3*omega*z)+c4*cos(4*omega*z)+c5*cos(5*omega*z)+c6*cos(6*omega*z)+c7*cos(7*omega*z)+c8*cos(8*omega*z)+c9*cos(9*omega*z)+c10*cos(10*omega*z)+...
    c11*cos(11*omega*z)+c12*cos(12*omega*z)+c13*cos(13*omega*z)+c14*cos(14*omega*z)+c15*cos(15*omega*z)+c16*cos(16*omega*z)+c17*cos(17*omega*z)+c18*cos(18*omega*z)+c19*cos(19*omega*z)+c20*cos(20*omega*z)+...
    c21*cos(21*omega*z)+c22*cos(22*omega*z)+c23*cos(23*omega*z)+c24*cos(24*omega*z)+c25*cos(25*omega*z)+c26*cos(26*omega*z)+c27*cos(27*omega*z)+c28*cos(28*omega*z)+c29*cos(29*omega*z)+c30*cos(30*omega*z)+...
    c31*cos(31*omega*z)+c32*cos(32*omega*z)+c33*cos(33*omega*z)+c34*cos(34*omega*z)+c35*cos(35*omega*z)+c36*cos(36*omega*z)+c37*cos(37*omega*z)+c38*cos(38*omega*z)+c39*cos(39*omega*z)+c40*cos(40*omega*z)+...
    d1*sin(omega*z)+d2*sin(2*omega*z)+d3*sin(3*omega*z)+d4*sin(4*omega*z)+d5*sin(5*omega*z)+d6*sin(6*omega*z)+d7*sin(7*omega*z)+d8*sin(8*omega*z)+d9*sin(9*omega*z)+d10*sin(10*omega*z)+...
    d11*sin(11*omega*z)+d12*sin(12*omega*z)+d13*sin(13*omega*z)+d14*sin(14*omega*z)+d15*sin(15*omega*z)+d16*sin(16*omega*z)+d17*sin(17*omega*z)+d18*sin(18*omega*z)+d19*sin(19*omega*z)+d20*sin(20*omega*z)+...
    d21*sin(21*omega*z)+d22*sin(22*omega*z)+d23*sin(23*omega*z)+d24*sin(24*omega*z)+d25*sin(25*omega*z)+d26*sin(26*omega*z)+d27*sin(27*omega*z)+d28*sin(28*omega*z)+d29*sin(29*omega*z)+d30*sin(30*omega*z)+...
    d31*sin(31*omega*z)+d32*sin(32*omega*z)+d33*sin(33*omega*z)+d34*sin(34*omega*z)+d35*sin(35*omega*z)+d36*sin(36*omega*z)+d37*sin(37*omega*z)+d38*sin(38*omega*z)+d39*sin(39*omega*z)+d40*sin(40*omega*z);


figure,plot(z,v,'gV',z,v,'-.',z,v6,'r+',z,v7,'k');grid;

figure,plot(z,((v-v7)./v)*100,'r-.');grid;
