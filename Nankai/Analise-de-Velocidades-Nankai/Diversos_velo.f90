subroutine velo_tti(modevel,nx,nz,dx,dz,tti,option)

real,intent(in),dimension(nz,nx)::modevel
integer,intent(in)::nx,nz,option
real,intent(in)::dz,dx
real,intent(inout),dimension(nz,nx)::tti

real::t0,profun,z1,x
integer::ix,iz,i,nz1

print*,modevel(1,1),dz/modevel(1,1)

do ix=1,nx

   x=dx+(ix-1)*dx
   t0=dz/modevel(1,ix)
   tti(1,ix)=t0
   z1=profun(x,option)
   nz1=int(z1/dz)

do iz=2,nz1
   do i=1,iz
tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))
   end do
tti(iz,ix)=(tti(iz,ix)+t0)*2.0
end do

do iz=nz1+1,nz
   do i=nz1,iz
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2
   end do
end do

end do

end subroutine velo_tti

subroutine velo_rms(modevel,tti,nx,nz,ntime,dx,dz,dt,velRMS,option)

real,intent(in),dimension(nz,nx)::modevel,tti
integer,intent(in)::nx,nz,ntime,option
real,intent(in)::dt,dx,dz
real,intent(inout),dimension(ntime,nx)::velRMS

integer::ix,iz,i,itime,it
real::x,z1,t

real,allocatable,dimension(:,:)::v_nmo,vN

allocate(v_nmo(nz,nx),vN(nz,nx))

do iz=1,nz
do ix=1,nx
v_nmo(iz,ix)=0.d0
vN(iz,ix)=0.d0
end do
end do

do ix=1,nx

   x=dx+(ix-1)*dx
   z1=profun(x,option)
   velRMS(1,ix)=modevel(1,ix)
   nz1=int(z1/dz)

do iz=2,nz1

   do i=1,iz
   if(i.eq.1) then
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)**2)*dt
   else
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)**2)
   end if
   end do

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix))

itime=int(tti(iz,ix)/dt)

if(itime<=ntime) then
velRMS(itime,ix)=vN(iz,ix)
else
end if

do it=2,itime
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime,ix)-velRMS(1,ix))*((it*dt-dt)/(itime*dt-dt))
end do

do it=itime+1,ntime
   velRMS(it,ix)=3500.0
end do

end do

end do

open(25,file='vel_rms_depth.dat',status='unknown') 
write(25,'(f10.5)') ((vN(j,i),j=1,nz),i=1,nx)
close(25,status='keep')

deallocate(v_nmo,vN)

end subroutine velo_rms

function profun(x,option) result(prof)

implicit none

real,intent(in)::x
integer,intent(in)::option
real::prof 

if(option.eq.1)then

if(x.gt.1305.and.x.lt.2405) then
prof=2000.d0+250.d0*(exp(-(x-1875.d0)**2/100000))
else
prof=2000
end if

else if(option.eq.2) then
prof=2000.d0
else if(option.eq.3) then
prof=2100.d0-0.176d0*x
else if(option.eq.4)then
prof=1000.d0      !600*(1+exp(-(x-600)**2/100000))+x/6

else if(option.eq.5) then
if(x.gt.1305.and.x.lt.2405) then
prof=2400.d0-250.d0*(exp(-(x-1875.d0)**2/100000))
else
prof=2400.d0
end if

else if(option.eq.6) then

    if(x.gt.2500) then
    prof=2000-0.5773*(x-2500);!%1501+x/2+(x/2-468)*max(0,(x-938));%a1=0.1245;a2=0.417/grau=-0.1763
    else if(x.le.1000) then
    prof=2000+0.5773*(x-1000)
    else if(x.gt.1000.and.x.le.2500) then
    prof=2000;
    end if;
else

end if

end function profun

function presal(x,option) result(prof)

implicit none

real,intent(in)::x
integer,intent(in)::option
real::prof 

if(option.eq.1)then        !!====== TOPO DO SAL ======!!

if(x.gt.5000.and.x.lt.17000.or.x.gt.19000.and.x.lt.25000) then
prof=4000.d0-1200.d0*(exp(-(x-8000.d0)**2/1000000))-1400.d0*(exp(-(x-22500.d0)**2/1000000))
else
prof=4000
end if

else if(option.eq.2) then  !!====== LÂMINA D'ÁGUA ======!!
prof=2000.d0

else if(option.eq.3) then  !!====== TOPO DO EMBASAMENTO ======!!

if(x.gt.0.and.x.lt.10000) then
prof=8000
else if(x.ge.10000.and.x.lt.20000) then
prof=9200.d0-0.12d0*x
else if (x.ge.20000) then
prof=6800
else
end if

else if(option.eq.4)then   !!====== BASE DO SAL ======!!
prof=6000.d0      

else if(option.eq.5) then
if(x.gt.700.and.x.lt.1780) then
prof=1400.d0-250.d0*(exp(-(x-1250.d0)**2/100000))
else
prof=1400.d0
end if

else if(option.eq.6) then

    if(x.gt.2500) then
    prof=1600-0.5773*(x-2500)
    else if(x.le.1000) then
    prof=1600+0.5773*(x-1000)
    else if(x.gt.1000.and.x.le.2500) then
    prof=1600;
    end if;
else

end if

end function presal

subroutine velo_tti_layers(modevel,nx,nz,dx,dz,tti,option)

real,intent(in),dimension(nz,nx)::modevel
integer,intent(in)::nx,nz,option
real,intent(in)::dz,dx
real,intent(inout),dimension(nz,nx)::tti

real::t0,presal,z1,z2,z3,z4,x
integer::ix,iz,i,nz1,nz2,nz3,nz4

print*,modevel(1,1),dz/modevel(1,1)

do ix=1,nx  ! Laço-x

   x=dx+(ix-1)*dx
   t0=dz/modevel(1,ix)
   tti(1,ix)=t0

   z1=presal(x,2) ! lamina d'agua
   z2=presal(x,1) ! topo do sal
   z3=presal(x,4) ! base do sal
   z4=presal(x,3) ! topo do embasamento

   nz1=int(z1/dz)
   nz2=int(z2/dz)
   nz3=int(z3/dz)
   nz4=int(z4/dz)

do iz=2,nz1
   do i=1,iz
tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))
   end do
tti(iz,ix)=(tti(iz,ix)+t0)*2.0
end do

do iz=nz1+1,nz2
   do i=nz1,iz
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0
   end do
   tti(iz,ix)=(tti(iz,ix)+tti(nz1,ix))*1.0
end do

do iz=nz2+1,nz3
   do i=nz2,iz
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0
   end do
   tti(iz,ix)=(tti(iz,ix)+tti(nz2,ix))
end do

do iz=nz3+1,nz4
   do i=nz3,iz
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0
   end do
   tti(iz,ix)=(tti(iz,ix)+tti(nz3,ix))
end do

do iz=nz4+1,nz
   do i=nz4,iz
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0
   end do
   tti(iz,ix)=(tti(iz,ix)+tti(nz4,ix))
end do

end do  ! Fim do laço-x

end subroutine velo_tti_layers

subroutine velo_rms_layers(modevel,tti,nx,nz,ntime,dx,dz,dt,velRMS,option)

real,intent(in),dimension(nz,nx)::modevel,tti
integer,intent(in)::nx,nz,ntime,option
real,intent(in)::dt,dx,dz
real,intent(inout),dimension(ntime,nx)::velRMS

integer::ix,iz,i,itime1,itime2,itime3,itime4,it,nz1,nz2,nz3,nz4
real::x,z1,z2,z3,z4,t,presal

real,allocatable,dimension(:,:)::v_nmo,vN

allocate(v_nmo(nz,nx),vN(nz,nx))

do iz=1,nz
do ix=1,nx
v_nmo(iz,ix)=0.d0
vN(iz,ix)=0.d0
end do
end do

do ix=1,nx  ! Laço-x

   x=dx+(ix-1)*dx

   z1=presal(x,2) ! lamina d'agua
   z2=presal(x,1) ! topo do sal
   z3=presal(x,4) ! base do sal
   z4=presal(x,3) ! topo do embasamento

   velRMS(1,ix)=modevel(1,ix)

   nz1=int(z1/dz)
   nz2=int(z2/dz)
   nz3=int(z3/dz)
   nz4=int(z4/dz)

do iz=2,nz1 ! Inicio camada-1

   do i=1,iz
   if(i.eq.1) then
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)**2)*dt
   else
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)**2)
   end if
   end do

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix))

itime1=int(tti(iz,ix)/dt)

if(itime1<=ntime) then
velRMS(itime1,ix)=vN(iz,ix)
else
end if

do it=2,itime1
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime1,ix)-velRMS(1,ix))*((it*dt-dt)/(itime1*dt-dt))
end do

end do ! Fim camada-1

do iz=nz1+1,nz2 ! Inicio camada-2

   do i=1,iz
   if(i.eq.1) then
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)**2)*dt
   else
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)**2)
   end if
   end do

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix))

itime2=int(tti(iz,ix)/dt)

if(itime2<=ntime) then
velRMS(itime2,ix)=vN(iz,ix)
else
end if

do it=itime1+1,itime2
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime2,ix)-velRMS(1,ix))*((it*dt-dt)/(itime2*dt-dt))
end do

end do ! Fim camada-2

do iz=nz2+1,nz3 ! Inicio camada-3

   do i=1,iz
   if(i.eq.1) then
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)**2)*dt
   else
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)**2)
   end if
   end do

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix))

itime3=int(tti(iz,ix)/dt)

if(itime3<=ntime) then
velRMS(itime3,ix)=vN(iz,ix)
else
end if

do it=itime2+1,itime3
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime3,ix)-velRMS(1,ix))*((it*dt-dt)/(itime3*dt-dt))
end do

end do ! Fim camada-3

do iz=nz3+1,nz4 ! Inicio camada-4

   do i=1,iz
   if(i.eq.1) then
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)**2)*dt
   else
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)**2)
   end if
   end do

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix))

itime4=int(tti(iz,ix)/dt)

if(itime4<=ntime) then
velRMS(itime4,ix)=vN(iz,ix)
else
end if

do it=itime3+1,itime4
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime4,ix)-velRMS(1,ix))*((it*dt-dt)/(itime4*dt-dt))
end do

end do ! Fim camada-4

do iz=nz4+1,nz ! inicio camada-5

do i=1,iz
   if(i.eq.1) then
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)**2)*dt
   else
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)**2)
   end if
   end do

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix))

itime5=int(tti(iz,ix)/dt)

if(itime5<=ntime) then
velRMS(itime5,ix)=vN(iz,ix)
else
end if

do it=itime4+1,ntime
   velRMS(it,ix)=velRMS(itime4,ix)+(velRMS(itime5,ix)-velRMS(itime4,ix))*((it*dt-itime4*dt)/(itime5*dt-itime4*dt))
end do

end do !Fim camada-5

end do  ! Fim do laço-x

open(25,file='vel_rms_depth.dat',status='unknown') 
write(25,'(f10.5)') ((vN(j,i),j=1,nz),i=1,nx)
close(25,status='keep')

deallocate(v_nmo,vN)

end subroutine velo_rms_layers
